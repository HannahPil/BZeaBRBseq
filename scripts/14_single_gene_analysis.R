# ==============================================================================
# INDIVIDUAL GENE EQTL PLOTS
# minimal setup for plotting expression of selected genes by introgression status
# required files in working directory:
#   results_list_new_name.rds
#   Zea_mays_counts.txt
#   metadata.csv
#   Zea_mays.gtf
# ==============================================================================

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# ##############################################################################
# ##                                                                          ##
# ##   >>> CHANGE THIS GENE ID TO ANALYZE A DIFFERENT GENE <<<                ##
# ##                                                                          ##
# ##############################################################################

FOCUS_GENE <- "Zm00001eb012750"

# partner genes for co-expression scatter plots against FOCUS_GENE
PARTNER_GENES <- c(
  "Zm00001eb399680",
  "Zm00001eb078420",
  "Zm00001eb179680"
)

# partner genes for "mixed-only" co-expression plots
PARTNER_GENES_MIXED <- c(
  "Zm00001eb399680",
  "Zm00001eb179680"
)

# ##############################################################################

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "individual_gene_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load data ------------------------------------

data_dir <- "data"

teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
counts   <- read.delim(file.path(data_dir, "Zea_mays_counts.txt"), check.names = FALSE, row.names = 1)
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

teogeno <- teogeno[!duplicated(names(teogeno))]

# ------------------------------ align samples ---------------------------------

sample_df <- metadata |>
  dplyr::filter(
    sample_id %in% colnames(counts),
    plate %in% c(1, 2, 3, 4)                # keep only sequenced plates
  ) |>
  dplyr::transmute(
    sample_id,
    taxa = factor(taxa),
    plate = factor(plate),
    genotype_teogeno_key = dplyr::if_else(
      genotype == "B73",
      "B73.B",
      paste0(genotype, ".B")
    )
  ) |>
  dplyr::filter(genotype_teogeno_key %in% names(teogeno)) |>
  dplyr::arrange(sample_id)

counts <- counts[, as.character(sample_df$sample_id), drop = FALSE]

# --------------------------- normalize expression -----------------------------

# NOTE: using raw library-size CPM (not TMM-normalized) here intentionally.
# Consistent with the eQTL scripts (13a, 13b).
# The PCA/DE scripts (10, 11, 12) use edgeR TMM CPM instead.
lib_size <- colSums(counts)
expr_mat <- as.matrix(counts)
expr_mat <- log2((t(t(expr_mat) / lib_size)) * 1e6 + 1)

# ---------------------------- gene coordinates --------------------------------

gtf <- import(file.path(data_dir, "Zea_mays.gtf"))
gtf_df <- as.data.frame(gtf)

gene_coords <- gtf_df |>
  dplyr::filter(!is.na(gene_id)) |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    chr = as.character(seqnames[1]),
    start = min(start),
    end = max(end),
    .groups = "drop"
  ) |>
  dplyr::filter(gene_id %in% rownames(expr_mat))

# ----------------------- introgression segments only --------------------------

seg_df <- purrr::imap_dfr(teogeno, ~{
  .x |>
    dplyr::filter(V4 == "Introgression") |>
    dplyr::transmute(
      key = .y,
      chr = V1,
      start = as.integer(V2),
      end = as.integer(V3)
    )
})

segs_gr <- GRanges(
  seqnames = seg_df$chr,
  ranges = IRanges(seg_df$start, seg_df$end),
  key = seg_df$key
)

key_to_samples <- split(sample_df$sample_id, sample_df$genotype_teogeno_key)

# ------------------------------ colors ----------------------------------------

taxa_colors <- c(
  "B73"  = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837",
  "Dura" = "#609bfe",
  "Nabo" = "#609bfe",
  "Mesa" = "#609bfe",
  "Chal" = "#609bfe",
  "Nobo" = "#609bfe"
)

# --------------------------- genes to plot ------------------------------------

genes_to_plot <- FOCUS_GENE

# --------------------------- plotting function --------------------------------

plot_one_gene <- function(target_gene) {
  
  if (!target_gene %in% gene_coords$gene_id) {
    stop(paste("gene not found in gene_coords:", target_gene))
  }
  
  if (!target_gene %in% rownames(expr_mat)) {
    stop(paste("gene not found in expression matrix:", target_gene))
  }
  
  gene_row <- gene_coords |>
    dplyr::filter(gene_id == target_gene)
  
  gene_gr <- GRanges(
    seqnames = gene_row$chr,
    ranges = IRanges(gene_row$start, gene_row$end),
    gene_id = gene_row$gene_id
  )
  
  hits <- findOverlaps(segs_gr, gene_gr)
  overlapping_keys <- unique(segs_gr$key[queryHits(hits)])
  
  gene_geno <- rep(0, nrow(sample_df))
  names(gene_geno) <- sample_df$sample_id
  
  for (k in overlapping_keys) {
    s_ids <- key_to_samples[[k]]
    if (!is.null(s_ids)) {
      gene_geno[as.character(s_ids)] <- 1
    }
  }
  
  gene_expr <- expr_mat[target_gene, sample_df$sample_id]
  
  plot_df <- tibble(
    sample_id = sample_df$sample_id,
    Genotype = factor(
      ifelse(gene_geno == 1, "Teosinte Introgression", "B73 Background"),
      levels = c("B73 Background", "Teosinte Introgression")
    ),
    Expression = as.numeric(gene_expr)
  ) |>
    dplyr::left_join(
      sample_df |>
        dplyr::select(sample_id, taxa),
      by = "sample_id"
    ) |>
    dplyr::mutate(
      point_color = ifelse(Genotype == "Teosinte Introgression", as.character(taxa), "B73")
    )
  
  group_stats <- plot_df |>
    dplyr::group_by(Genotype) |>
    dplyr::summarise(
      mean_expr = mean(Expression),
      sd_expr = sd(Expression),
      .groups = "drop"
    )
  
  p_gene <- ggplot(plot_df, aes(x = Genotype, y = Expression)) +
    geom_jitter(
      aes(fill = point_color),
      shape = 21,
      size = 3,
      alpha = 0.75,
      color = "grey30",
      position = position_jitter(width = 0.2)
    ) +
    geom_errorbar(
      data = group_stats,
      aes(
        x = Genotype,
        ymin = mean_expr - sd_expr,
        ymax = mean_expr + sd_expr
      ),
      width = 0,
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = group_stats,
      aes(x = Genotype, y = mean_expr),
      shape = 16,
      size = 4,
      color = "black",
      inherit.aes = FALSE
    ) +
    scale_x_discrete(labels = c(
      "B73 Background" = "B73",
      "Teosinte Introgression" = "Teo"
    )) +
    scale_fill_manual(
      values = taxa_colors,
      name = "teosinte taxa"
    ) +
    labs(
      x = "genotype",
      y = "normalized expression (log2 cpm)",
      title = paste("Effect of introgression on", target_gene),
      subtitle = "Teo points colored by taxa; black = mean ± SD"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  print(p_gene)
  
  ggsave(
    file.path(out_dir, paste0("eQTL_", target_gene, "_taxa.png")),
    p_gene,
    width = 5,
    height = 5
  )
  
  return(plot_df)
}

# ------------------------------ run plots -------------------------------------

plot_results <- purrr::map(genes_to_plot, plot_one_gene)
names(plot_results) <- genes_to_plot

# ##############################################################################
# ##  SISTER-LINE ANALYSIS                                                    ##
# ##  Compare expression of FOCUS_GENE in teo-carriers vs their BC1 sisters   ##
# ##  (related lines sharing the same BC1 parent that lack the introgression) ##
# ##############################################################################

# ---- metadata with pedigree info (exclude checks) ----------------------------

meta_ped <- metadata |>
  dplyr::filter(
    sample_id %in% colnames(expr_mat),
    plate %in% c(1, 2, 3, 4),
    !genotype %in% c("B73", "Purple Check")
  ) |>
  dplyr::mutate(
    bc1_family = substr(genotype, 1, 13)
  )

# ---- introgression status for FOCUS_GENE per sample --------------------------

focus_row <- gene_coords |> dplyr::filter(gene_id == FOCUS_GENE)
focus_gr  <- GRanges(
  seqnames = focus_row$chr,
  ranges   = IRanges(focus_row$start, focus_row$end)
)

hits_focus <- findOverlaps(segs_gr, focus_gr)
overlapping_keys <- unique(segs_gr$key[queryHits(hits_focus)])

focus_geno <- rep(0L, nrow(sample_df))
names(focus_geno) <- sample_df$sample_id
for (k in overlapping_keys) {
  s_ids <- key_to_samples[[k]]
  if (!is.null(s_ids)) focus_geno[as.character(s_ids)] <- 1L
}

meta_ped <- meta_ped |>
  dplyr::mutate(
    has_teo = focus_geno[sample_id],
    allele  = ifelse(has_teo == 1, "Teosinte", "B73")
  )

# ---- list all lines carrying teosinte for FOCUS_GENE -------------------------

teo_carriers_df <- meta_ped |>
  dplyr::filter(has_teo == 1) |>
  dplyr::select(sample_id, genotype, taxa, bc1_family) |>
  dplyr::distinct()

cat("\n=== Lines carrying teosinte introgression at", FOCUS_GENE, "===\n")
print(as.data.frame(teo_carriers_df), row.names = FALSE)

write.csv(
  teo_carriers_df,
  file.path(out_dir, paste0(FOCUS_GENE, "_teo_carriers.csv")),
  row.names = FALSE
)

# ---- find sister lines (same BC1 family, no introgression) -------------------

teo_families <- unique(teo_carriers_df$bc1_family)

sister_df <- meta_ped |>
  dplyr::filter(
    bc1_family %in% teo_families,   # same BC1 parent
    has_teo == 0                     # but lacks introgression at FOCUS_GENE
  ) |>
  dplyr::select(sample_id, genotype, taxa, bc1_family) |>
  dplyr::distinct()

cat("\n=== Sister lines (same BC1, B73 allele at", FOCUS_GENE, ") ===\n")
print(as.data.frame(sister_df), row.names = FALSE)
cat("\nFamilies with teo + sister pairs:", length(intersect(teo_families, sister_df$bc1_family)), "\n")

# ---- comprehensive carrier + sister table with expression --------------------

carrier_sister_df <- meta_ped |>
  dplyr::filter(bc1_family %in% teo_families) |>
  dplyr::mutate(
    expression = as.numeric(expr_mat[FOCUS_GENE, sample_id]),
    role       = ifelse(has_teo == 1, "Teo", "Sis")
  ) |>
  dplyr::select(
    bc1_family, role, sample_id, genotype, taxa, expression
  ) |>
  dplyr::arrange(bc1_family, dplyr::desc(role), dplyr::desc(expression))

write.csv(
  carrier_sister_df,
  file.path(out_dir, paste0(FOCUS_GENE, "_carriers_and_sisters.csv")),
  row.names = FALSE
)

cat("\n=== Carriers and sisters with expression (sorted by family) ===\n")
print(as.data.frame(carrier_sister_df), row.names = FALSE)

# ---- build paired expression table ------------------------------------------

paired_df <- meta_ped |>
  dplyr::filter(bc1_family %in% teo_families) |>
  dplyr::mutate(
    expression = as.numeric(expr_mat[FOCUS_GENE, sample_id])
  )

# ---- effect size plot: teo vs B73 within BC1 families ------------------------

# family means for effect size
family_means <- paired_df |>
  dplyr::group_by(bc1_family, allele) |>
  dplyr::summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from  = allele,
    values_from = c(mean_expr, n)
  ) |>
  dplyr::filter(!is.na(mean_expr_Teosinte), !is.na(mean_expr_B73)) |>
  dplyr::mutate(
    effect = mean_expr_Teosinte - mean_expr_B73
  ) |>
  dplyr::left_join(
    meta_ped |> dplyr::select(bc1_family, taxa) |> dplyr::distinct(),
    by = "bc1_family"
  )

cat("\n=== Effect sizes (teo - B73) within BC1 families ===\n")
print(as.data.frame(family_means |> dplyr::select(bc1_family, taxa, effect, n_Teosinte, n_B73)), row.names = FALSE)

write.csv(
  family_means,
  file.path(out_dir, paste0(FOCUS_GENE, "_sister_effect_sizes.csv")),
  row.names = FALSE
)

# ---- waterfall bar chart (effect per family, sorted) -------------------------

# drop families with NA taxa (no teosinte annotation)
family_means_sorted <- family_means |>
  dplyr::filter(!is.na(taxa)) |>
  dplyr::arrange(effect) |>
  dplyr::mutate(bc1_family = factor(bc1_family, levels = bc1_family))

p_waterfall <- ggplot(family_means_sorted, aes(x = bc1_family, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = taxa), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = taxa_colors, name = "Taxa") +
  labs(
    x = "BC1 family",
    y = "Effect size (Teo - B73, log2 CPM)",
    title = paste("Introgression effect by family:", FOCUS_GENE),
    subtitle = "Positive = teosinte allele upregulates"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title  = element_text(face = "bold")
  )

print(p_waterfall)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_waterfall.png")),
  p_waterfall,
  width = max(6, nrow(family_means_sorted) * 0.4 + 2),
  height = 5
)

# ---- faceted jitter plots per family (mean ± SD style) -----------------------
# Same style as the main B73-vs-Teo jitter plot: filled circles colored by taxa,
# black mean dot, and SD error bar per allele group within each BC1 family.

# helper function to build the faceted jitter plot
build_facet_jitter <- function(df_fam, fam_levels, show_legend = TRUE,
                               title_suffix = "", nrow = NULL) {

  df_fam <- df_fam |>
    dplyr::mutate(
      bc1_family = factor(bc1_family, levels = fam_levels),
      allele     = factor(allele, levels = c("B73", "Teosinte"))
    )

  fstats <- df_fam |>
    dplyr::group_by(bc1_family, allele) |>
    dplyr::summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      sd_expr   = sd(expression, na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    ) |>
    dplyr::mutate(sd_expr = ifelse(is.na(sd_expr), 0, sd_expr))

  p <- ggplot(df_fam, aes(x = allele, y = expression)) +
    geom_jitter(
      aes(fill = taxa),
      shape    = 21,
      size     = 3,
      alpha    = 0.75,
      color    = "grey30",
      position = position_jitter(width = 0.2)
    ) +
    geom_errorbar(
      data = fstats,
      aes(
        x    = allele,
        ymin = mean_expr - sd_expr,
        ymax = mean_expr + sd_expr
      ),
      width     = 0,
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = fstats,
      aes(x = allele, y = mean_expr),
      shape = 16,
      size  = 4,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ bc1_family, nrow = nrow) +
    scale_fill_manual(values = taxa_colors, name = "teosinte taxa") +
    scale_x_discrete(labels = c("B73" = "B73", "Teosinte" = "Teo")) +
    labs(
      x = "genotype",
      y = "normalized expression (log2 cpm)",
      title = paste0("Per-family expression: ", FOCUS_GENE, title_suffix),
      subtitle = "Each panel = one BC1 family; black = mean ± SD"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold"),
      strip.text  = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold")
    )

  if (!show_legend) p <- p + theme(legend.position = "none")

  return(p)
}

# sort families by effect size; keep only non-NA taxa
family_order <- levels(family_means_sorted$bc1_family)

paired_df_fam <- paired_df |>
  dplyr::filter(bc1_family %in% teo_families, !is.na(taxa)) |>
  dplyr::filter(bc1_family %in% family_means_sorted$bc1_family)

# -- all families plot --
p_facet <- build_facet_jitter(paired_df_fam, family_order)
print(p_facet)

n_fam <- length(unique(paired_df_fam$bc1_family))
ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_faceted_families.png")),
  p_facet,
  width  = min(16, max(6, ceiling(sqrt(n_fam)) * 3)),
  height = min(14, max(4, ceiling(n_fam / ceiling(sqrt(n_fam))) * 3))
)

# -- top 5 families by absolute effect size (no legend) --
top5_families <- family_means_sorted |>
  dplyr::arrange(dplyr::desc(abs(effect))) |>
  dplyr::slice_head(n = 5) |>
  dplyr::arrange(effect) |>
  dplyr::pull(bc1_family) |>
  as.character()

paired_df_top5 <- paired_df_fam |>
  dplyr::filter(bc1_family %in% top5_families)

p_top5 <- build_facet_jitter(
  paired_df_top5,
  top5_families,
  show_legend  = FALSE,
  title_suffix = " (top 5)",
  nrow         = 1
)
print(p_top5)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_faceted_top5.png")),
  p_top5,
  width  = 16,
  height = 4
)

# ---------------------- pairwise co-expression plots --------------------------

target_gene <- FOCUS_GENE
partner_genes <- PARTNER_GENES

get_gene_intro_status <- function(target_gene) {
  
  if (!target_gene %in% gene_coords$gene_id) {
    stop(paste("gene not found in gene_coords:", target_gene))
  }
  
  gene_row <- gene_coords |>
    dplyr::filter(gene_id == target_gene)
  
  gene_gr <- GRanges(
    seqnames = gene_row$chr,
    ranges = IRanges(gene_row$start, gene_row$end),
    gene_id = gene_row$gene_id
  )
  
  hits <- findOverlaps(segs_gr, gene_gr)
  overlapping_keys <- unique(segs_gr$key[queryHits(hits)])
  
  gene_geno <- rep(0, nrow(sample_df))
  names(gene_geno) <- sample_df$sample_id
  
  for (k in overlapping_keys) {
    s_ids <- key_to_samples[[k]]
    if (!is.null(s_ids)) {
      gene_geno[as.character(s_ids)] <- 1
    }
  }
  
  tibble(
    sample_id = sample_df$sample_id,
    intro_status = gene_geno
  )
}

plot_gene_pair <- function(y_gene, x_gene) {
  
  if (!y_gene %in% rownames(expr_mat)) {
    stop(paste("y gene not found in expression matrix:", y_gene))
  }
  
  if (!x_gene %in% rownames(expr_mat)) {
    stop(paste("x gene not found in expression matrix:", x_gene))
  }
  
  y_status <- get_gene_intro_status(y_gene) |>
    dplyr::rename(y_intro = intro_status)
  
  x_status <- get_gene_intro_status(x_gene) |>
    dplyr::rename(x_intro = intro_status)
  
  plot_df <- tibble(
    sample_id = sample_df$sample_id,
    x_expr = as.numeric(expr_mat[x_gene, sample_df$sample_id]),
    y_expr = as.numeric(expr_mat[y_gene, sample_df$sample_id])
  ) |>
    dplyr::left_join(y_status, by = "sample_id") |>
    dplyr::left_join(x_status, by = "sample_id") |>
    dplyr::mutate(
      pair_status = dplyr::case_when(
        x_intro == 0 & y_intro == 0 ~ "B73 / B73",
        x_intro == 1 & y_intro == 1 ~ "Teo / Teo",
        TRUE ~ "Mixed"
      ),
      pair_status = factor(
        pair_status,
        levels = c("B73 / B73", "Mixed", "Teo / Teo")
      )
    )
  
  p_pair <- ggplot(plot_df, aes(x = x_expr, y = y_expr, color = pair_status)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      x = paste0(x_gene, " expression (log2 cpm)"),
      y = paste0(y_gene, " expression (log2 cpm)"),
      title = paste("Co-expression:", y_gene, "vs", x_gene),
      subtitle = "point color = introgression status for the two genes"
    ) +
    scale_color_manual(
      values = c(
        "B73 / B73" = "gray50",
        "Mixed" = "goldenrod3",
        "Teo / Teo" = "forestgreen"
      )
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  print(p_pair)
  
  ggsave(
    file.path(out_dir, paste0("coexpression_", y_gene, "_vs_", x_gene, ".png")),
    p_pair,
    width = 6,
    height = 5
  )
  
  return(plot_df)
}

pair_plot_results <- purrr::map(partner_genes, ~plot_gene_pair(target_gene, .x))
names(pair_plot_results) <- partner_genes

# ---------------- mixed-only co-expression plots by which gene is teo ---------

target_gene_mixed <- FOCUS_GENE
partner_genes_mixed <- PARTNER_GENES_MIXED

plot_gene_pair_mixed_only <- function(y_gene, x_gene) {
  
  if (!y_gene %in% rownames(expr_mat)) {
    stop(paste("y gene not found in expression matrix:", y_gene))
  }
  
  if (!x_gene %in% rownames(expr_mat)) {
    stop(paste("x gene not found in expression matrix:", x_gene))
  }
  
  y_status <- get_gene_intro_status(y_gene) |>
    dplyr::rename(y_intro = intro_status)
  
  x_status <- get_gene_intro_status(x_gene) |>
    dplyr::rename(x_intro = intro_status)
  
  plot_df <- tibble(
    sample_id = sample_df$sample_id,
    x_expr = as.numeric(expr_mat[x_gene, sample_df$sample_id]),
    y_expr = as.numeric(expr_mat[y_gene, sample_df$sample_id])
  ) |>
    dplyr::left_join(y_status, by = "sample_id") |>
    dplyr::left_join(x_status, by = "sample_id") |>
    dplyr::mutate(
      mixed_group = dplyr::case_when(
        y_intro == 1 & x_intro == 0 ~ paste0(y_gene, " teo"),
        y_intro == 0 & x_intro == 1 ~ paste0(x_gene, " teo"),
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(mixed_group)) |>
    dplyr::mutate(
      mixed_group = factor(
        mixed_group,
        levels = c(paste0(y_gene, " teo"), paste0(x_gene, " teo"))
      )
    )
  
  p_pair <- ggplot(plot_df, aes(x = x_expr, y = y_expr, color = mixed_group)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      x = paste0(x_gene, " expression (log2 cpm)"),
      y = paste0(y_gene, " expression (log2 cpm)"),
      title = paste("Mixed samples only:", y_gene, "vs", x_gene),
      subtitle = "blue = y-axis gene is teo, purple = x-axis gene is teo"
    ) +
    scale_color_manual(
      values = setNames(
        c("dodgerblue3", "purple3"),
        c(paste0(y_gene, " teo"), paste0(x_gene, " teo"))
      )
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  print(p_pair)
  
  ggsave(
    file.path(out_dir, paste0("coexpression_mixed_only_", y_gene, "_vs_", x_gene, ".png")),
    p_pair,
    width = 6,
    height = 5
  )
  
  return(plot_df)
}

pair_plot_results_mixed_only <- purrr::map(
  partner_genes_mixed,
  ~plot_gene_pair_mixed_only(target_gene_mixed, .x)
)

names(pair_plot_results_mixed_only) <- partner_genes_mixed

# ##############################################################################
# ##                                                                          ##
# ##   NDH COMPLEX VOLCANO PLOT                                               ##
# ##   Effect of teosinte allele on expression of all NDH complex genes       ##
# ##                                                                          ##
# ##############################################################################

library(ggrepel)

ndh_file <- file.path(data_dir, "NDH_complex_genes_B73v5_COMPLETE.tsv")
ndh_genes <- read.delim(ndh_file, stringsAsFactors = FALSE)

# keep only nuclear-encoded genes with valid IDs present in expression matrix
ndh_nuclear <- ndh_genes |>
  dplyr::filter(Genome == "nuclear", Gene_ID != "-", Gene_ID %in% rownames(expr_mat))

cat(paste("NDH nuclear genes in expression data:", nrow(ndh_nuclear), "\n"))

# ---- compute per-gene effect of teosinte allele -----------------------------

# For each gene, determine introgression status per sample, then compare
# expression between teo-carriers and B73-carriers at that locus.

compute_teo_effect <- function(gene_id, gene_coords_df, segs_gr,
                               key_to_samples, sample_df, expr_mat) {

  if (!gene_id %in% gene_coords_df$gene_id) return(NULL)

  g <- gene_coords_df |> dplyr::filter(gene_id == !!gene_id)
  g_gr <- GRanges(seqnames = g$chr, ranges = IRanges(g$start, g$end))

  hits <- findOverlaps(segs_gr, g_gr)
  overlapping_keys <- unique(segs_gr$key[queryHits(hits)])

  geno_vec <- rep(0L, nrow(sample_df))
  names(geno_vec) <- sample_df$sample_id
  for (k in overlapping_keys) {
    s_ids <- key_to_samples[[k]]
    if (!is.null(s_ids)) geno_vec[as.character(s_ids)] <- 1L
  }

  expr_vals <- expr_mat[gene_id, names(geno_vec)]
  teo_expr <- expr_vals[geno_vec == 1]
  b73_expr <- expr_vals[geno_vec == 0]

  n_teo <- length(teo_expr)
  n_b73 <- length(b73_expr)

  if (n_teo < 2 || n_b73 < 2) return(NULL)

  tt <- t.test(teo_expr, b73_expr)

  data.frame(
    gene_id      = gene_id,
    mean_teo     = mean(teo_expr),
    mean_b73     = mean(b73_expr),
    log2FC       = mean(teo_expr) - mean(b73_expr),  # already log2 CPM
    pvalue       = tt$p.value,
    n_teo        = n_teo,
    n_b73        = n_b73,
    stringsAsFactors = FALSE
  )
}

ndh_effects <- purrr::map_dfr(
  ndh_nuclear$Gene_ID,
  ~compute_teo_effect(.x, gene_coords, segs_gr, key_to_samples, sample_df, expr_mat)
)

# merge subunit annotations back
ndh_effects <- ndh_effects |>
  dplyr::left_join(
    ndh_nuclear |> dplyr::select(Gene_ID, Subcomplex, Subunit),
    by = c("gene_id" = "Gene_ID")
  ) |>
  dplyr::mutate(
    neg_log10p   = -log10(pvalue),
    padj         = p.adjust(pvalue, method = "BH"),
    neg_log10padj = -log10(padj),
    is_focus     = gene_id == FOCUS_GENE,
    sig          = padj < 0.05,
    label        = Subunit
  )

cat("\n=== NDH complex teosinte allele effects ===\n")
print(
  as.data.frame(
    ndh_effects |>
      dplyr::select(Subcomplex, Subunit, gene_id, log2FC, pvalue, padj, n_teo, n_b73) |>
      dplyr::arrange(pvalue)
  ),
  row.names = FALSE
)

write.csv(
  ndh_effects,
  file.path(out_dir, "NDH_complex_teo_effects.csv"),
  row.names = FALSE
)

# ---- volcano plot ------------------------------------------------------------

subcomplex_colors <- c(
  "SubA"       = "#E64B35",
  "SubB"       = "#4DBBD5",
  "SubE"       = "#00A087",
  "SubL"       = "#3C5488",
  "SubM"       = "#F39B7F",
  "PSI-linker" = "#8491B4"
)

p_volcano <- ggplot(ndh_effects, aes(x = log2FC, y = neg_log10padj)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(
    aes(fill = Subcomplex, size = is_focus),
    shape = 21,
    color = "grey30",
    alpha = 0.8
  ) +
  scale_fill_manual(values = subcomplex_colors, name = "Subcomplex") +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
  geom_text_repel(
    aes(label = label),
    size          = 3.5,
    fontface      = "bold",
    max.overlaps  = 30,
    segment.color = "grey60",
    segment.size  = 0.3,
    box.padding   = 0.4,
    point.padding = 0.3
  ) +
  labs(
    x = "Effect of teosinte allele (log2 FC)",
    y = expression(-log[10]~adjusted~italic(p)-value),
    title = "NDH complex: teosinte introgression effect on expression",
    subtitle = paste0(
      sum(ndh_effects$sig), " / ", nrow(ndh_effects),
      " genes significant (BH-adjusted p < 0.05)  |  ",
      "Focus gene (PnsL1/PPL2) highlighted"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )

print(p_volcano)

ggsave(
  file.path(out_dir, "NDH_complex_volcano.png"),
  p_volcano,
  width = 9,
  height = 7
)

cat("Saved: NDH_complex_volcano.png\n")

# ---- jitter plots for significant NDH genes ----------------------------------

ndh_sig_genes <- ndh_effects |>
  dplyr::filter(sig) |>
  dplyr::pull(gene_id)

cat(paste("\nGenerating jitter plots for", length(ndh_sig_genes), "significant NDH genes...\n"))

ndh_jitter_results <- purrr::map(ndh_sig_genes, plot_one_gene)
names(ndh_jitter_results) <- ndh_sig_genes

cat("Saved jitter plots for all significant NDH genes.\n")