#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 1b (LOCAL) — R_div / R_cons for fbxl1 (Rubén memo §6.3)
#
# Reads window counts produced by scripts/FBX_analysis1_count_windows.R on HPC
# (now written directly into data/FBX_*). Classifies each sample by fbxl1
# (Zm00001eb375600) introgression status using the same logic as A4. Computes:
#   R_cons = mean(Teo) / mean(B73) in sense_terminalCDS_conserved  (bias-free)
#   R_div  = mean(Teo) / mean(B73) in sense_3UTR_divergent          (bias-exposed)
#   R_div / R_cons  = mapping-bias diagnostic (H1)
#
# Interpretation (from memo):
#   R_div / R_cons ~ 1              -> no bias, expression difference is real (H1 rejected)
#   R_div / R_cons < 1              -> reads lost specifically in divergent window
#                                       ratio itself is a direct estimate of bias
#   effect gone once bias removed    -> H1 alone (no cis-regulatory effect)
#   effect persists after bias fix   -> H2 or H3 candidate; re-estimate from CDS
#
# Outputs:
#   output/FBX_analysis1_ratio/window_ratios_summary.csv
#   output/FBX_analysis1_ratio/per_sample_window_cpm.png
#   output/FBX_analysis1_ratio/ratio_by_window.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

FOCUS_GENE <- "Zm00001eb012750"    # actually PPL2 — wait, this is fbxl1
FOCUS_GENE <- "Zm00001eb375600"    # fbxl1

data_dir <- "data"
out_dir  <- file.path("output", "FBX_analysis1_ratio")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. HPC window counts + library sizes --------------------------------
win_norm <- read.csv(file.path(data_dir, "FBX_window_counts_normalized.csv"),
                     check.names = FALSE, row.names = 1)
win_raw  <- read.csv(file.path(data_dir, "FBX_window_counts.csv"),
                     check.names = FALSE, row.names = 1)
lib_size <- read.csv(file.path(data_dir, "FBX_library_sizes.csv"),
                     stringsAsFactors = FALSE)

cat("Windows counted:\n")
print(rownames(win_norm))
cat("\nSamples:", ncol(win_norm), "\n")

# ---- 2. Classify samples by fbxl1 introgression status -------------------
# Mirrors SG_single_gene_analysis.R sample_df logic (keeps B73 controls)
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

sample_df <- metadata |>
  dplyr::filter(sample_id %in% colnames(win_norm),
                plate %in% c(1, 2, 3, 4)) |>
  dplyr::transmute(
    sample_id, genotype, taxa = as.character(taxa), plate,
    genotype_teogeno_key = dplyr::if_else(genotype == "B73",
                                          "B73.B", paste0(genotype, ".B"))
  ) |>
  dplyr::filter(genotype_teogeno_key %in% names(teogeno) |
                  genotype == "B73") |>
  dplyr::arrange(sample_id)

gtf <- import(file.path(data_dir, "external", "Zea_mays.gtf"))
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

cat("\nfbxl1 introgression classification:\n")
print(table(sample_df$Genotype))

# ---- 3. Long-form data: per sample x window ------------------------------
common <- intersect(colnames(win_norm), sample_df$sample_id)
if (length(common) < nrow(sample_df)) {
  warning("Some samples missing from window counts: ",
          length(sample_df$sample_id) - length(common))
}
win_norm <- win_norm[, common, drop = FALSE]
win_raw  <- win_raw[,  common, drop = FALSE]
sample_df <- sample_df |> dplyr::filter(sample_id %in% common)

long_df <- as.data.frame(win_norm) |>
  tibble::rownames_to_column("window") |>
  tidyr::pivot_longer(-window, names_to = "sample_id", values_to = "cpm") |>
  dplyr::left_join(
    as.data.frame(win_raw) |>
      tibble::rownames_to_column("window") |>
      tidyr::pivot_longer(-window, names_to = "sample_id", values_to = "raw"),
    by = c("window", "sample_id")
  ) |>
  dplyr::inner_join(sample_df |>
                      dplyr::select(sample_id, Genotype, taxa),
                    by = "sample_id")

# ---- 4. R_cons, R_div, R_div/R_cons --------------------------------------
summarise_window <- function(df, window_name) {
  sub <- df |> dplyr::filter(window == window_name)
  b73 <- sub |> dplyr::filter(Genotype == "B73") |> dplyr::pull(cpm)
  teo <- sub |> dplyr::filter(Genotype == "Teo") |> dplyr::pull(cpm)
  data.frame(
    window       = window_name,
    n_B73        = length(b73),
    n_Teo        = length(teo),
    mean_cpm_B73 = mean(b73),
    mean_cpm_Teo = mean(teo),
    median_cpm_B73 = median(b73),
    median_cpm_Teo = median(teo),
    ratio_Teo_over_B73 = mean(teo) / mean(b73)
  )
}

summary_tbl <- purrr::map_dfr(rownames(win_norm), ~summarise_window(long_df, .x))
cat("\nPer-window summary (Teo vs B73 at fbxl1):\n")
print(summary_tbl, row.names = FALSE, digits = 4)

# extract R_cons, R_div, and the diagnostic
R_cons <- summary_tbl$ratio_Teo_over_B73[
  summary_tbl$window == "sense_terminalCDS_conserved"]
R_div  <- summary_tbl$ratio_Teo_over_B73[
  summary_tbl$window == "sense_3UTR_divergent"]
R_ratio <- R_div / R_cons

cat("\n=========================================================================\n")
cat("MEMO §6.3 KEY NUMBERS\n")
cat("=========================================================================\n")
cat(sprintf("  R_cons  (Teo/B73 in conserved terminal CDS window) = %.3f\n", R_cons))
cat(sprintf("  R_div   (Teo/B73 in divergent 3'UTR window)        = %.3f\n", R_div))
cat(sprintf("  R_div / R_cons                                     = %.3f\n", R_ratio))
cat("\nInterpretation:\n")
if (abs(R_ratio - 1) < 0.15) {
  cat("  R_div/R_cons is close to 1 -> read recovery is equally good in\n",
      "  both windows. Expression difference is REAL abundance, not\n",
      "  mapping bias. H1 REJECTED. Move to H2/H3.\n", sep = "")
} else if (R_ratio < 0.85) {
  cat("  R_div/R_cons is < 1 -> reads are being lost specifically in the\n",
      "  divergent 3'UTR window. That drop is a direct estimate of the\n",
      "  mapping bias. The conserved-window ratio (R_cons = ",
      round(R_cons, 3), ") is the bias-corrected biological effect.\n",
      "  If R_cons is close to 1, the fbxl1 result is ENTIRELY mapping bias.\n",
      "  If R_cons is well below 1, bias INFLATED a real effect.\n", sep = "")
} else if (R_ratio > 1.15) {
  cat("  R_div/R_cons > 1 -> unexpected; teo reads recover BETTER in\n",
      "  divergent window than in conserved window. Something odd is happening;\n",
      "  investigate before drawing conclusions.\n", sep = "")
}

# save summary
write.csv(summary_tbl, file.path(out_dir, "window_ratios_summary.csv"),
          row.names = FALSE)
cat("\nSummary written to: ", file.path(out_dir, "window_ratios_summary.csv"),
    "\n", sep = "")

# ---- 5. Plots ------------------------------------------------------------
# Per-sample CPM in each of the three sense windows, split by genotype
sense_windows <- c("sense_terminalCDS_conserved", "sense_3UTR_divergent",
                   "antisense_test_intron4")
plot_df <- long_df |>
  dplyr::filter(window %in% sense_windows) |>
  dplyr::mutate(window = factor(window, levels = sense_windows))

p1 <- ggplot(plot_df, aes(x = Genotype, y = cpm + 0.1, fill = Genotype)) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.5) +
  scale_y_log10() +
  facet_wrap(~ window, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("B73" = "#03bec4", "Teo" = "#d62728"),
                    guide = "none") +
  labs(
    title = "fbxl1 windows — per-sample CPM by introgression status",
    subtitle = "(y-axis log10; +0.1 pseudocount for the log)",
    x = NULL,
    y = "reads per million (window)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "per_sample_window_cpm.png"),
       p1, width = 10, height = 5, dpi = 200)

# Ratio bar chart — one bar per window
ratio_df <- summary_tbl |>
  dplyr::filter(window %in% sense_windows) |>
  dplyr::mutate(
    window = factor(window, levels = sense_windows),
    label = sprintf("%.2f", ratio_Teo_over_B73)
  )

p2 <- ggplot(ratio_df, aes(x = window, y = ratio_Teo_over_B73, fill = window)) +
  geom_col(width = 0.6, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = label), vjust = -0.3, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("sense_terminalCDS_conserved" = "#4a90d9",
                               "sense_3UTR_divergent"        = "#d94848",
                               "antisense_test_intron4"      = "#888888"),
                    guide = "none") +
  labs(
    title = "Teo / B73 mean CPM ratio per window",
    subtitle = sprintf("R_div / R_cons = %.3f (memo diagnostic)", R_ratio),
    x = NULL, y = "mean(Teo) / mean(B73)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out_dir, "ratio_by_window.png"),
       p2, width = 8, height = 5, dpi = 200)

cat("\nPlots written to:\n")
cat("  ", file.path(out_dir, "per_sample_window_cpm.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "ratio_by_window.png"), "\n", sep = "")
