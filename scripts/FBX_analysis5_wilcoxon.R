#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 5 (LOCAL) — per-taxon Wilcoxon test on CDS-window CPM
# for Rubén reply memo §6.3 / §10 item 5.
#
# Uses the featureCounts-based CPM values (sense_terminalCDS_conserved window
# in FBX_window_counts_normalized.csv). Rank-sum test is essentially a median
# comparison; Rubén specifically asked for it stratified by donor taxon
# because running across all 39 carriers pools taxa with likely different
# read-quality profiles.
#
# For each Teo donor taxon:
#   - Wilcoxon rank-sum vs the B73 background samples on CDS-window CPM
#   - Median and mean CPM, ratio to B73 median/mean
# Plus overall Teo (all carriers pooled) vs B73.
#
# Outputs:
#   output/FBX_analysis5_wilcoxon/per_taxon_wilcoxon.csv
#   output/FBX_analysis5_wilcoxon/wilcoxon_forest.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

FOCUS_GENE   <- "Zm00001eb375600"    # fbxl1
CDS_WINDOW   <- "sense_terminalCDS_conserved"

data_dir <- "data"
out_dir  <- file.path("output", "FBX_analysis5_wilcoxon")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Per-sample CDS-window CPM ---------------------------------------
win_norm <- read.csv(file.path(data_dir, "FBX_window_counts_normalized.csv"),
                     check.names = FALSE, row.names = 1)
stopifnot(CDS_WINDOW %in% rownames(win_norm))
cds_cpm <- unlist(win_norm[CDS_WINDOW, ])
cat("Samples with CDS-window CPM: ", length(cds_cpm), "\n", sep = "")

# ---- 2. Classify by fbxl1 introgression (same logic as analyses 1/3) ----
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

sample_df <- metadata |>
  dplyr::filter(sample_id %in% names(cds_cpm), plate %in% c(1, 2, 3, 4)) |>
  dplyr::transmute(
    sample_id, genotype, taxa = as.character(taxa), plate,
    genotype_teogeno_key = dplyr::if_else(genotype == "B73",
                                          "B73.B", paste0(genotype, ".B"))
  ) |>
  dplyr::filter(genotype_teogeno_key %in% names(teogeno) | genotype == "B73") |>
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
gene_geno <- setNames(rep(0L, nrow(sample_df)), sample_df$sample_id)
for (k in overlap_keys) {
  s <- key_to_samples[[k]]
  if (!is.null(s)) gene_geno[as.character(s)] <- 1L
}
sample_df$Genotype <- factor(
  ifelse(gene_geno == 1, "Teo", "B73"), levels = c("B73", "Teo"))
sample_df$cds_cpm <- cds_cpm[sample_df$sample_id]

cat("\nGroup sizes:\n"); print(table(sample_df$Genotype))

# Lab-shorthand -> donor taxon (from accession-code prefix in metadata.csv)
taxa_to_donor <- c(
  Bals = "parviglumis",     Chal = "mexicana",
  Dura = "mexicana",        Hueh = "huehuetenangensis",
  Mesa = "mexicana",        Nobo = "mexicana",
  Zdip = "diploperennis",   Zlux = "luxurians"
)

# ---- 3. Wilcoxon per taxon ----------------------------------------------
b73_vec <- sample_df$cds_cpm[sample_df$Genotype == "B73"]

run_wilcox <- function(name, teo_vec) {
  if (length(teo_vec) < 2) {
    return(data.frame(comparison = name, n_teo = length(teo_vec),
                      n_b73 = length(b73_vec),
                      median_teo = median(teo_vec),
                      median_b73 = median(b73_vec),
                      ratio_median = median(teo_vec) / median(b73_vec),
                      mean_teo = mean(teo_vec),
                      mean_b73 = mean(b73_vec),
                      ratio_mean = mean(teo_vec) / mean(b73_vec),
                      W = NA, p_value = NA))
  }
  w <- suppressWarnings(wilcox.test(teo_vec, b73_vec, alternative = "two.sided"))
  data.frame(
    comparison   = name,
    n_teo        = length(teo_vec),
    n_b73        = length(b73_vec),
    median_teo   = median(teo_vec),
    median_b73   = median(b73_vec),
    ratio_median = median(teo_vec) / median(b73_vec),
    mean_teo     = mean(teo_vec),
    mean_b73     = mean(b73_vec),
    ratio_mean   = mean(teo_vec) / mean(b73_vec),
    W            = unname(w$statistic),
    p_value      = w$p.value
  )
}

teo_all <- sample_df$cds_cpm[sample_df$Genotype == "Teo"]
overall <- run_wilcox("All Teo carriers (pooled)", teo_all)

per_taxon <- purrr::map_dfr(
  sort(unique(sample_df$taxa[sample_df$Genotype == "Teo"])),
  function(tx) {
    vec <- sample_df$cds_cpm[sample_df$Genotype == "Teo" & sample_df$taxa == tx]
    df  <- run_wilcox(tx, vec)
    df$donor <- unname(taxa_to_donor[tx])
    df
  }
)
overall$donor <- NA_character_
result <- dplyr::bind_rows(overall, per_taxon) |>
  dplyr::mutate(p_bonferroni = pmin(1, p_value * nrow(per_taxon)))

cat("\nPer-taxon Wilcoxon rank-sum vs B73 background:\n")
print(result, digits = 4, row.names = FALSE)
write.csv(result, file.path(out_dir, "per_taxon_wilcoxon.csv"), row.names = FALSE)

# ---- 4. Forest-style plot -----------------------------------------------
plot_df <- result |>
  dplyr::filter(comparison != "All Teo carriers (pooled)") |>
  dplyr::mutate(
    signif = dplyr::case_when(
      is.na(p_value)          ~ "n<2",
      p_bonferroni < 0.05     ~ "signif (Bonferroni)",
      p_value < 0.05          ~ "nominal p<0.05",
      TRUE                    ~ "n.s."
    ),
    label_p = ifelse(is.na(p_value), "n<2",
                     ifelse(p_value < 0.001, "p<0.001",
                            sprintf("p=%.3f", p_value)))
  )

col_sig <- c(
  "signif (Bonferroni)" = "#d62728",
  "nominal p<0.05"      = "#f4a261",
  "n.s."                = "#03bec4",
  "n<2"                 = "grey70"
)

p_forest <- ggplot(plot_df,
                   aes(x = ratio_median,
                       y = reorder(comparison, ratio_median),
                       colour = signif)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(label_p, "  (n=", n_teo, ")")),
            hjust = -0.15, size = 3.2, colour = "grey20", show.legend = FALSE) +
  scale_colour_manual(values = col_sig, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title    = "Per-taxon Wilcoxon on CDS-window CPM (Teo vs B73 background)",
    subtitle = sprintf("Overall pooled: median ratio = %.2f, p = %.3g (n_Teo=%d, n_B73=%d)",
                       overall$ratio_median, overall$p_value,
                       overall$n_teo, overall$n_b73),
    x = "median(Teo) / median(B73) in CDS window",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(out_dir, "wilcoxon_forest.png"),
       p_forest, width = 9, height = 5, dpi = 200)

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "per_taxon_wilcoxon.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "wilcoxon_forest.png"), "\n", sep = "")
