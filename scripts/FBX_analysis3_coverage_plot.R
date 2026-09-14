#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 3 (LOCAL) — per-base coverage profile of fbxl1 terminal exon
# Rubén reply memo §5 (decisive test) and §6.1 (taxon-stratified CDS ratio).
#
# INPUTS (copy from HPC to local data/FBX_analysis3/):
#   data/FBX_analysis3/depth_matrix.tsv         (from analysis3_coverage/)
#   data/FBX_analysis1/library_sizes.csv        (already local from analysis 1)
#
# TWO OUTPUTS:
#   output/FBX_analyses/analysis3_coverage/coverage_profile.png
#       Per-base mean library-normalised depth (carriers vs B73 background)
#       across chr9:17933600-17934260. Distinguishes:
#         (i)  lower teo expression        -> same shape, scaled down
#         (ii) shorter teo 3'UTR (geometry)-> shape shifted upstream
#         (iii) mapping loss in divergent  -> drop only in the 3'UTR window
#
#   output/FBX_analyses/analysis3_coverage/cds_ratio_by_taxon.{png,csv}
#       Rubén §6.1 falsification test for the geometry model. If the 3.08x
#       CDS-window ratio tracks donor 3'UTR length across carrier taxa,
#       geometry explains it. If it is uniform, real expression difference.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

FOCUS_GENE   <- "Zm00001eb375600"    # fbxl1
POLY_A_POS   <- 17934180L            # B73 poly(A) site (end of 3'UTR window)
CDS_START    <- 17933693L; CDS_END   <- 17933902L   # 210 bp
UTR_START    <- 17933903L; UTR_END   <- 17934180L   # 278 bp
REGION_START <- 17933600L; REGION_END <- 17934260L  # samtools depth region

data_dir <- "data"
depth_dir <- if (file.exists(file.path(data_dir, "FBX_analysis3", "depth_matrix.tsv"))) {
  file.path(data_dir, "FBX_analysis3")
} else {
  data_dir
}
lib_dir <- if (file.exists(file.path(data_dir, "FBX_analysis1", "library_sizes.csv"))) {
  file.path(data_dir, "FBX_analysis1")
} else {
  data_dir
}
out_dir <- file.path("output", "FBX_analyses", "analysis3_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Depth matrix from:", depth_dir, "\n")
cat("Library sizes from:", lib_dir, "\n")

# ---- 1. Depth matrix + library sizes -------------------------------------
depth <- read.table(file.path(depth_dir, "depth_matrix.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE)
stopifnot(all(c("chr", "pos") %in% colnames(depth)))
cat("Depth matrix: ", nrow(depth), " positions x ", ncol(depth) - 2,
    " samples\n", sep = "")

lib_size <- read.csv(file.path(lib_dir, "library_sizes.csv"),
                     stringsAsFactors = FALSE)
lib_vec <- setNames(lib_size$lib_size, lib_size$sample_id)

# ---- 2. Classify each sample by fbxl1 introgression (same as analysis 1) -
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

sample_df <- metadata |>
  dplyr::filter(sample_id %in% colnames(depth),
                plate %in% c(1, 2, 3, 4)) |>
  dplyr::transmute(
    sample_id, genotype, taxa = as.character(taxa), plate,
    genotype_teogeno_key = dplyr::if_else(genotype == "B73",
                                          "B73.B", paste0(genotype, ".B"))
  ) |>
  dplyr::filter(genotype_teogeno_key %in% names(teogeno) |
                  genotype == "B73") |>
  dplyr::arrange(sample_id)

gtf <- import(file.path(data_dir, "reference", "Zea_mays.gtf"))
gene_coords <- as.data.frame(gtf) |>
  dplyr::filter(!is.na(gene_id)) |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(chr = as.character(seqnames[1]),
                   start = min(start), end = max(end), .groups = "drop")
seg_df <- purrr::imap_dfr(teogeno, ~{
  .x |> dplyr::filter(V4 == "Introgression") |>
    dplyr::transmute(key = .y, chr = V1,
                     start = as.integer(V2), end = as.integer(V3))
})
segs_gr <- GRanges(seqnames = seg_df$chr,
                   ranges = IRanges(seg_df$start, seg_df$end),
                   key = seg_df$key)
gr <- gene_coords |> dplyr::filter(gene_id == FOCUS_GENE)
focus_gr <- GRanges(seqnames = gr$chr,
                    ranges = IRanges(gr$start, gr$end))
overlap_keys <- unique(segs_gr$key[queryHits(findOverlaps(segs_gr, focus_gr))])
key_to_samples <- split(sample_df$sample_id, sample_df$genotype_teogeno_key)

gene_geno <- rep(0L, nrow(sample_df))
names(gene_geno) <- sample_df$sample_id
for (k in overlap_keys) {
  s <- key_to_samples[[k]]
  if (!is.null(s)) gene_geno[as.character(s)] <- 1L
}
sample_df$Genotype <- factor(
  ifelse(gene_geno == 1, "Teo", "B73"),
  levels = c("B73", "Teo")
)
cat("\nfbxl1 introgression classification:\n"); print(table(sample_df$Genotype))
cat("\nTeo carriers by taxa:\n")
print(sample_df |> dplyr::filter(Genotype == "Teo") |>
        dplyr::count(taxa, sort = TRUE))

# ---- 3. Long depth table with normalisation & Genotype -------------------
common <- intersect(colnames(depth)[-(1:2)], sample_df$sample_id)
common <- intersect(common, names(lib_vec))
if (length(common) < nrow(sample_df)) {
  warning("Dropping samples missing from either depth or lib_size: ",
          nrow(sample_df) - length(common))
}
sample_df <- sample_df |> dplyr::filter(sample_id %in% common)

# per-sample per-position normalised depth (depth per million mapped)
dp <- depth[, c("chr", "pos", common), drop = FALSE]
dp_long <- dp |>
  tidyr::pivot_longer(-c(chr, pos), names_to = "sample_id",
                      values_to = "depth") |>
  dplyr::inner_join(sample_df |> dplyr::select(sample_id, Genotype, taxa),
                    by = "sample_id") |>
  dplyr::mutate(depth_norm = depth / lib_vec[sample_id] * 1e6,
                pos_rel    = pos - POLY_A_POS)

# ---- 4. Coverage profile (Rubén reply memo Fig 1c reproduction) ----------
cov_by_group <- dp_long |>
  dplyr::group_by(Genotype, pos, pos_rel) |>
  dplyr::summarise(
    mean_depth_norm   = mean(depth_norm),
    median_depth_norm = median(depth_norm),
    n = dplyr::n(),
    .groups = "drop"
  )

# also normalise each group to its own max, so shape vs area is separable
cov_by_group <- cov_by_group |>
  dplyr::group_by(Genotype) |>
  dplyr::mutate(rel_depth = mean_depth_norm / max(mean_depth_norm)) |>
  dplyr::ungroup()

col_geno <- c("B73" = "#03bec4", "Teo" = "#d62728")

# absolute (mean normalised) coverage
p_abs <- ggplot(cov_by_group,
                aes(x = pos_rel, y = mean_depth_norm, colour = Genotype)) +
  annotate("rect", xmin = CDS_START - POLY_A_POS, xmax = CDS_END - POLY_A_POS,
           ymin = -Inf, ymax = Inf, fill = "#4a90d9", alpha = 0.10) +
  annotate("rect", xmin = UTR_START - POLY_A_POS, xmax = UTR_END - POLY_A_POS,
           ymin = -Inf, ymax = Inf, fill = "#d94848", alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = col_geno) +
  labs(
    title    = "fbxl1 terminal exon — mean per-base depth (library-normalised)",
    subtitle = "Blue band: CDS window (conserved). Red band: 3'UTR window (divergent). Dashed line: B73 poly(A) site.",
    x = "bp from B73 poly(A) site",
    y = "mean depth per million mapped reads"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# shape-only (each group scaled to its own max) -> tells geometry from expression
p_shape <- ggplot(cov_by_group,
                  aes(x = pos_rel, y = rel_depth, colour = Genotype)) +
  annotate("rect", xmin = CDS_START - POLY_A_POS, xmax = CDS_END - POLY_A_POS,
           ymin = -Inf, ymax = Inf, fill = "#4a90d9", alpha = 0.10) +
  annotate("rect", xmin = UTR_START - POLY_A_POS, xmax = UTR_END - POLY_A_POS,
           ymin = -Inf, ymax = Inf, fill = "#d94848", alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = col_geno) +
  labs(
    title    = "Shape only — each group rescaled to its own maximum",
    subtitle = "Same shape scaled -> pure expression diff. Shape shifted upstream -> shorter teo UTR (geometry). Drop only in red -> mapping loss.",
    x = "bp from B73 poly(A) site",
    y = "relative depth (group max = 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# stack the two panels
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  p_combined <- p_abs / p_shape
} else {
  p_combined <- p_abs
}

ggsave(file.path(out_dir, "coverage_profile.png"),
       p_combined, width = 10, height = 8, dpi = 200)
write.csv(cov_by_group,
          file.path(out_dir, "coverage_profile_data.csv"),
          row.names = FALSE)

# ---- 5. Rubén §6.1 — taxon-stratified CDS-window ratio -------------------
# Compute per-sample CDS-window CPM (mean depth over 17933693-17933902,
# normalised by library size then averaged), then group by donor taxa.

cds_per_sample <- dp_long |>
  dplyr::filter(pos >= CDS_START, pos <= CDS_END) |>
  dplyr::group_by(sample_id, Genotype, taxa) |>
  dplyr::summarise(cds_mean_depth_norm = mean(depth_norm), .groups = "drop")

utr_per_sample <- dp_long |>
  dplyr::filter(pos >= UTR_START, pos <= UTR_END) |>
  dplyr::group_by(sample_id, Genotype, taxa) |>
  dplyr::summarise(utr_mean_depth_norm = mean(depth_norm), .groups = "drop")

b73_cds_mean <- mean(cds_per_sample$cds_mean_depth_norm[
  cds_per_sample$Genotype == "B73"])
b73_utr_mean <- mean(utr_per_sample$utr_mean_depth_norm[
  utr_per_sample$Genotype == "B73"])

# per-taxon summary among Teo carriers
teo_taxon <- cds_per_sample |>
  dplyr::filter(Genotype == "Teo") |>
  dplyr::group_by(taxa) |>
  dplyr::summarise(
    n_carriers          = dplyr::n(),
    mean_cds            = mean(cds_mean_depth_norm),
    median_cds          = median(cds_mean_depth_norm),
    ratio_mean_over_B73_cds   = mean(cds_mean_depth_norm) / b73_cds_mean,
    ratio_median_over_B73_cds = median(cds_mean_depth_norm) / b73_cds_mean,
    .groups = "drop"
  ) |>
  dplyr::left_join(
    utr_per_sample |>
      dplyr::filter(Genotype == "Teo") |>
      dplyr::group_by(taxa) |>
      dplyr::summarise(
        mean_utr = mean(utr_mean_depth_norm),
        ratio_mean_over_B73_utr = mean(utr_mean_depth_norm) / b73_utr_mean,
        .groups = "drop"),
    by = "taxa"
  ) |>
  dplyr::arrange(dplyr::desc(ratio_mean_over_B73_cds))

# Rubén's 3'UTR length table from GFF3s (mexicana + parviglumis only;
# diploperennis/luxurians/huehuetenangensis: NA — flag those explicitly).
utr_lengths_bp <- tribble(
  ~taxa_key,                    ~utr_bp,
  "parviglumis",                215L,        # Zv-TIL01
  "mexicana",                   285L,        # Zx-TIL18 (Zx-TIL25 = 293)
  "diploperennis",              NA_integer_, # not verified
  "luxurians",                  NA_integer_,
  "huehuetenangensis",          NA_integer_
)
teo_taxon <- teo_taxon |>
  dplyr::mutate(taxa_lower = tolower(taxa)) |>
  dplyr::mutate(donor_utr_bp = utr_lengths_bp$utr_bp[
    match(sapply(taxa_lower, function(x)
      utr_lengths_bp$taxa_key[which(sapply(utr_lengths_bp$taxa_key,
                                           function(k) grepl(k, x)))[1]]),
      utr_lengths_bp$taxa_key)])

cat("\nPer-donor-taxon CDS-window summary (Rubén §6.1 test):\n")
print(teo_taxon, n = Inf)
write.csv(teo_taxon, file.path(out_dir, "cds_ratio_by_taxon.csv"),
          row.names = FALSE)

# bar chart: CDS ratio per taxon, with observed 3.08x reference and
# geometry-model prediction range (2.12 - 4.06x) shaded.
p_taxon <- ggplot(teo_taxon,
                  aes(x = reorder(taxa, ratio_mean_over_B73_cds),
                      y = ratio_mean_over_B73_cds)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.12, ymax = 4.06,
           fill = "#ffcc99", alpha = 0.35) +
  geom_hline(yintercept = 1,    linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 3.08, linetype = "dotted", colour = "black") +
  geom_col(aes(fill = !is.na(donor_utr_bp)), width = 0.6, colour = "black") +
  geom_text(aes(label = paste0("n=", n_carriers)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(`TRUE` = "#03bec4", `FALSE` = "grey70"),
                    labels = c(`TRUE` = "3'UTR length known",
                               `FALSE` = "not verified"),
                    name   = NULL) +
  labs(
    title = "CDS-window Teo/B73 ratio, stratified by donor taxon",
    subtitle = "Orange band: geometry-model prediction 2.12-4.06x. Dotted: pooled 3.08x. Dashed: 1.",
    x = NULL, y = "mean(Teo CDS) / mean(B73 CDS)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "cds_ratio_by_taxon.png"),
       p_taxon, width = 8, height = 5, dpi = 200)

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "coverage_profile.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "coverage_profile_data.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "cds_ratio_by_taxon.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "cds_ratio_by_taxon.csv"), "\n", sep = "")
