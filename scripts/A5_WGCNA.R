# ==============================================================================
# WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA)
# Identifies co-expression modules, correlates them with introgression status
# and taxa, and finds hub genes in the focus gene's module.
# required files in data/:
#   results_list_new_name.rds
#   Zea_mays_counts.txt
#   metadata.csv
#   Zea_mays.gtf
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# allow multi-threading for WGCNA (uses all available cores)
allowWGCNAThreads()

# ##############################################################################
# ##                                                                          ##
# ##   >>> CONFIGURATION <<<                                                  ##
# ##                                                                          ##
# ##   FOCUS_GENE : gene to highlight in module analysis                      ##
# ##   N_TOP_GENES: number of most-variable genes to include (default 5000)   ##
# ##   SOFT_POWER : set to NULL for automatic detection, or override (e.g. 6) ##
# ##                                                                          ##
# ##############################################################################

FOCUS_GENE  <- "Zm00001eb012750"
N_TOP_GENES <- 5000
SOFT_POWER  <- NULL        # NULL = auto-detect; set to integer to override

# ##############################################################################

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "WGCNA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 1: DATA LOADING
# ==============================================================================

data_dir <- "data"

teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
counts   <- read.delim(file.path(data_dir, "Zea_mays_counts.txt"), check.names = FALSE, row.names = 1)
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

teogeno <- teogeno[!duplicated(names(teogeno))]

# ---- recode taxa to lowercase subspecies names -------------------------------
# The 5 mexicana subpopulations (Dura/Nabo/Mesa/Chal/Nobo) collapse to a
# single "mexicana" group.
taxa_recode_map <- c(
  "Bals" = "parviglumis",
  "Zdip" = "diploperennis",
  "Hueh" = "huehuetenanguensis",
  "Zlux" = "luxurians",
  "Dura" = "mexicana",
  "Nabo" = "mexicana",
  "Mesa" = "mexicana",
  "Chal" = "mexicana",
  "Nobo" = "mexicana",
  "B73"  = "B73"
)
recode_taxa <- function(x) {
  mapped <- unname(taxa_recode_map[as.character(x)])
  ifelse(is.na(mapped), as.character(x), mapped)
}
metadata$taxa <- recode_taxa(metadata$taxa)

# ------------------------------ align samples ---------------------------------

sample_df <- metadata |>
  dplyr::filter(
    sample_id %in% colnames(counts),
    plate %in% c(1, 2, 3, 4)                # keep only sequenced plates
  ) |>
  dplyr::transmute(
    sample_id,
    genotype,
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
# Consistent with the eQTL scripts (13a, 13b, 14).
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
  "B73"                = "#03bec4",
  "parviglumis"        = "#f364e2",
  "diploperennis"      = "#f8756d",
  "huehuetenanguensis" = "#b69d00",
  "luxurians"          = "#00b837",
  "mexicana"           = "#609bfe"
)

# ==============================================================================
# SECTION 2: GENE FILTERING — keep top variable genes + FOCUS_GENE
# ==============================================================================

gene_var <- apply(expr_mat, 1, var)

# also require genes to be expressed in at least some samples
gene_expressed <- rowSums(expr_mat > 0) >= 10

# rank by variance among expressed genes
gene_var_filtered <- gene_var[gene_expressed]
var_threshold <- sort(gene_var_filtered, decreasing = TRUE)[min(N_TOP_GENES, length(gene_var_filtered))]

keep_genes <- names(gene_var_filtered[gene_var_filtered >= var_threshold])

# always include FOCUS_GENE even if below variance cutoff
if (!FOCUS_GENE %in% keep_genes && FOCUS_GENE %in% rownames(expr_mat)) {
  keep_genes <- c(keep_genes, FOCUS_GENE)
  cat("Note: FOCUS_GENE added despite low variance\n")
}

cat(paste("Genes kept for WGCNA:", length(keep_genes), "\n"))

# WGCNA wants samples in rows, genes in columns
datExpr <- t(expr_mat[keep_genes, ])

# check for problematic genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("Removed", sum(!gsg$goodSamples), "samples and",
      sum(!gsg$goodGenes), "genes with too many missing values\n")
}

cat(paste("Final matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n"))

# ==============================================================================
# SECTION 3: SOFT-THRESHOLD SELECTION
# ==============================================================================

# WGCNA overrides base::cor() with its own version. We must use WGCNA's cor
# during network construction, then restore base R's cor for downstream code.
cor <- WGCNA::cor

powers <- c(1:10, seq(12, 20, 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose     = 5,
  networkType = "unsigned"
)

# auto-select: first power where scale-free R^2 > 0.85
if (is.null(SOFT_POWER)) {
  r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  idx <- which(r2 > 0.85)[1]
  if (is.na(idx)) {
    SOFT_POWER <- 6
    cat("Warning: no power reached R^2 > 0.85; defaulting to power = 6\n")
  } else {
    SOFT_POWER <- powers[idx]
  }
}
cat(paste("Selected soft-threshold power:", SOFT_POWER, "\n"))

# ---- plot soft-threshold diagnostics -----------------------------------------

png(file.path(out_dir, "WGCNA_soft_threshold.png"), width = 10, height = 5, units = "in", res = 150)
par(mfrow = c(1, 2))

# scale-free topology fit
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (signed R^2)",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = 0.9,
  col = "red"
)
abline(h = 0.85, col = "red", lty = 2)

# mean connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = 0.9,
  col = "red"
)

dev.off()
cat("Saved: WGCNA_soft_threshold.png\n")

# ==============================================================================
# SECTION 4: NETWORK CONSTRUCTION & MODULE DETECTION
# ==============================================================================

net <- blockwiseModules(
  datExpr,
  power                = SOFT_POWER,
  TOMType              = "unsigned",
  minModuleSize        = 30,
  reassignThreshold    = 0,
  mergeCutHeight       = 0.25,
  numericLabels        = TRUE,
  pamRespectsDendro    = FALSE,
  saveTOMs             = FALSE,
  verbose              = 3
)

# convert numeric labels to color labels
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)

# restore base R cor() now that WGCNA network construction is done
cor <- stats::cor

cat(paste("Number of modules detected:", length(unique(moduleColors)) - 1, "(+ grey)\n"))
cat("Module sizes:\n")
print(sort(table(moduleColors), decreasing = TRUE))

# ==============================================================================
# SECTION 5: MODULE VISUALIZATION
# ==============================================================================

# ---- dendrogram with module colors ------------------------------------------

png(file.path(out_dir, "WGCNA_dendrogram.png"), width = 12, height = 6, units = "in", res = 150)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()
cat("Saved: WGCNA_dendrogram.png\n")

# ---- module eigengenes -------------------------------------------------------

MEs <- net$MEs
colnames(MEs) <- gsub("^ME", "", colnames(MEs))

# reorder by eigengene similarity
ME_diss <- 1 - cor(MEs)
ME_tree <- hclust(as.dist(ME_diss), method = "average")

png(file.path(out_dir, "WGCNA_eigengene_clustering.png"), width = 8, height = 5, units = "in", res = 150)
plot(ME_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()
cat("Saved: WGCNA_eigengene_clustering.png\n")

# ==============================================================================
# SECTION 6: MODULE-TRAIT CORRELATION
# ==============================================================================

# ---- build trait matrix ------------------------------------------------------

# introgression status at FOCUS_GENE (same logic as script 14)
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

# trait dataframe aligned to datExpr rows
trait_df <- sample_df |>
  dplyr::filter(sample_id %in% rownames(datExpr)) |>
  dplyr::mutate(
    allele = focus_geno[sample_id]      # 1 = teo, 0 = B73
  )

# one-hot encode taxa (drop first level to avoid collinearity)
taxa_onehot <- model.matrix(~ taxa - 1, data = trait_df)
colnames(taxa_onehot) <- gsub("^taxa", "taxa_", colnames(taxa_onehot))

# one-hot encode plate
plate_onehot <- model.matrix(~ plate - 1, data = trait_df)
colnames(plate_onehot) <- gsub("^plate", "plate_", colnames(plate_onehot))

trait_mat <- cbind(
  allele = trait_df$allele,
  taxa_onehot,
  plate_onehot
)
rownames(trait_mat) <- trait_df$sample_id

# align to datExpr row order
trait_mat <- trait_mat[rownames(datExpr), , drop = FALSE]

# ---- correlate module eigengenes with traits ---------------------------------

MEs_ordered <- MEs[rownames(datExpr), , drop = FALSE]

module_trait_cor <- cor(MEs_ordered, trait_mat, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# ---- heatmap -----------------------------------------------------------------

# text matrix for display (r + stars)
textMatrix <- paste0(
  signif(module_trait_cor, 2),
  ifelse(module_trait_pval < 0.001, "***",
    ifelse(module_trait_pval < 0.01, "**",
      ifelse(module_trait_pval < 0.05, "*", "")))
)
dim(textMatrix) <- dim(module_trait_cor)

png(file.path(out_dir, "WGCNA_module_trait_heatmap.png"),
    width = max(8, ncol(trait_mat) * 0.6 + 2),
    height = max(6, ncol(MEs) * 0.4 + 2),
    units = "in", res = 150)

par(mar = c(8, 10, 3, 2))
labeledHeatmap(
  Matrix    = module_trait_cor,
  xLabels   = colnames(trait_mat),
  yLabels   = colnames(MEs_ordered),
  ySymbols  = colnames(MEs_ordered),
  colorLabels = FALSE,
  colors    = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text  = 0.5,
  zlim      = c(-1, 1),
  main      = "Module-trait relationships"
)

dev.off()
cat("Saved: WGCNA_module_trait_heatmap.png\n")

# ==============================================================================
# SECTION 7: FOCUS GENE MODULE ANALYSIS
# ==============================================================================

focus_module <- moduleColors[FOCUS_GENE]
cat(paste("\n>>> FOCUS_GENE", FOCUS_GENE, "is in module:", focus_module, "<<<\n"))

# ---- module membership (kME) for all genes -----------------------------------

kME <- as.data.frame(cor(datExpr, MEs_ordered, use = "p"))
colnames(kME) <- paste0("kME_", colnames(kME))

# kME for the focus module specifically
focus_me_col <- paste0("kME_", which(labels2colors(0:(ncol(MEs) - 1)) == focus_module) - 1)

# safer approach: use the module eigengene directly
focus_ME <- MEs_ordered[, labels2colors(as.integer(colnames(MEs_ordered))) == focus_module, drop = FALSE]
if (ncol(focus_ME) == 0) {
  # fallback: find the column by matching color
  me_colors <- labels2colors(as.integer(colnames(MEs_ordered)))
  focus_col_idx <- which(me_colors == focus_module)
  if (length(focus_col_idx) > 0) {
    focus_ME <- MEs_ordered[, focus_col_idx, drop = FALSE]
  }
}

kME_focus <- cor(datExpr, focus_ME, use = "p")[, 1]
names(kME_focus) <- colnames(datExpr)

# ---- hub genes (top 20 by kME in the focus module) ---------------------------

module_genes <- names(moduleColors[moduleColors == focus_module])
kME_module <- kME_focus[module_genes]

hub_genes <- data.frame(
  gene_id = names(sort(kME_module, decreasing = TRUE)),
  kME     = sort(kME_module, decreasing = TRUE),
  is_focus_gene = names(sort(kME_module, decreasing = TRUE)) == FOCUS_GENE,
  stringsAsFactors = FALSE
) |>
  dplyr::slice_head(n = 20)

cat(paste("\nTop 20 hub genes in the", focus_module, "module:\n"))
print(hub_genes, row.names = FALSE)

write.csv(
  hub_genes,
  file.path(out_dir, paste0("WGCNA_", FOCUS_GENE, "_hub_genes.csv")),
  row.names = FALSE
)

# ---- full module assignments -------------------------------------------------

module_assignments <- data.frame(
  gene_id = names(moduleColors),
  module  = moduleColors,
  kME_own_module = kME_focus[names(moduleColors)],
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(module, dplyr::desc(kME_own_module))

write.csv(
  module_assignments,
  file.path(out_dir, "WGCNA_module_assignments.csv"),
  row.names = FALSE
)
cat("Saved: WGCNA_module_assignments.csv\n")

# ==============================================================================
# SECTION 8: GENE SIGNIFICANCE vs MODULE MEMBERSHIP PLOT
# ==============================================================================

# gene significance = correlation of each gene with allele (introgression status)
allele_vec <- trait_mat[, "allele"]
gene_significance <- cor(datExpr, allele_vec, use = "p")[, 1]

# build plot dataframe for genes in the focus module
gs_mm_df <- data.frame(
  gene_id = module_genes,
  kME     = kME_focus[module_genes],
  GS      = gene_significance[module_genes],
  is_focus = module_genes == FOCUS_GENE,
  stringsAsFactors = FALSE
)

# correlation between GS and MM
gs_mm_cor <- cor(gs_mm_df$kME, gs_mm_df$GS, use = "p")
gs_mm_pval <- cor.test(gs_mm_df$kME, gs_mm_df$GS)$p.value

p_gs_mm <- ggplot(gs_mm_df, aes(x = kME, y = GS)) +
  geom_point(
    aes(size = is_focus, color = is_focus),
    alpha = 0.6
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey50", "TRUE" = "red"),
    guide  = "none"
  ) +
  scale_size_manual(
    values = c("FALSE" = 1.5, "TRUE" = 4),
    guide  = "none"
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 0.8) +
  ggrepel::geom_text_repel(
    data = gs_mm_df |> dplyr::filter(is_focus),
    aes(label = gene_id),
    color = "red",
    fontface = "bold",
    size = 4,
    nudge_y = 0.05
  ) +
  labs(
    x = paste("Module Membership (kME) in", focus_module, "module"),
    y = paste("Gene Significance for allele at", FOCUS_GENE),
    title = paste("GS vs MM:", focus_module, "module"),
    subtitle = paste0("r = ", round(gs_mm_cor, 3),
                      ", p = ", signif(gs_mm_pval, 3),
                      " | ", length(module_genes), " genes")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p_gs_mm)

ggsave(
  file.path(out_dir, paste0("WGCNA_", FOCUS_GENE, "_GS_vs_MM.png")),
  p_gs_mm,
  width = 7,
  height = 6
)

# ==============================================================================
# SECTION 9: CROSS-REFERENCE WITH EQTL RESULTS (if available)
# ==============================================================================

# check if trans-eQTL results exist
trans_file <- file.path("output", "trans_by_source_chr", "BZea_trans_ALL_sources_combined.tsv")

if (file.exists(trans_file)) {
  trans_results <- read.delim(trans_file)

  # find genes in the focus module that are also trans-eQTL targets
  module_eqtl <- trans_results |>
    dplyr::filter(gene %in% module_genes | SNP %in% module_genes)

  if (nrow(module_eqtl) > 0) {
    cat(paste("\n>>>", nrow(module_eqtl),
              "trans-eQTL hits involve genes in the", focus_module, "module <<<\n"))

    write.csv(
      module_eqtl,
      file.path(out_dir, paste0("WGCNA_", FOCUS_GENE, "_module_eQTL_overlap.csv")),
      row.names = FALSE
    )
    cat("Saved: module-eQTL overlap table\n")
  } else {
    cat("No trans-eQTL hits overlap the focus module.\n")
  }
} else {
  cat("Trans-eQTL results not found; skipping cross-reference.\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("=======================================================\n")
cat("  WGCNA ANALYSIS COMPLETE\n")
cat("=======================================================\n")
cat(paste("  Focus gene:     ", FOCUS_GENE, "\n"))
cat(paste("  Module:          ", focus_module, "\n"))
cat(paste("  Genes in module: ", length(module_genes), "\n"))
cat(paste("  Total modules:   ", length(unique(moduleColors)) - 1, "(+ grey)\n"))
cat(paste("  Soft power:      ", SOFT_POWER, "\n"))
cat(paste("  Genes analyzed:  ", ncol(datExpr), "\n"))
cat(paste("  Samples:         ", nrow(datExpr), "\n"))
cat("=======================================================\n")
cat("Output files in:", out_dir, "\n")
