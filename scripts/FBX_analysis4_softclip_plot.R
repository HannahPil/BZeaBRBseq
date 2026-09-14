#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 4 (LOCAL) — soft-clip / start-position audit
# Rubén reply memo §6.2 falsification test.
#
# Reads the per-read TSV produced on HPC (writes into data/FBX_softclip_reads.tsv
# via git pull) and, per sample, computes:
#   - fraction of reads with any soft-clip
#   - mean trailing soft-clip length (softR) — the 3'-edge signature
#   - histogram of read start positions across the CDS window
#
# Splits samples by fbxl1 introgression (Teo carrier vs B73 background).
# Interpretation:
#   Real CDS reads:      start positions spread across the window;
#                        soft-clip rate similar to B73.
#   Displaced UTR reads: start positions cluster against the 3' edge;
#                        higher soft-clip rate (esp. softR) in Teo.
#
# Outputs:
#   output/FBX_analysis4_softclip/softclip_rate_by_group.png
#   output/FBX_analysis4_softclip/start_position_density.png
#   output/FBX_analysis4_softclip/softR_length_density.png
#   output/FBX_analysis4_softclip/per_sample_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

FOCUS_GENE <- "Zm00001eb375600"    # fbxl1
CDS_START  <- 17933693L
CDS_END    <- 17933902L

data_dir <- "data"
out_dir  <- file.path("output", "FBX_analysis4_softclip")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load per-read table ---------------------------------------------
reads <- read.table(file.path(data_dir, "FBX_softclip_reads.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Reads: ", nrow(reads), " across ",
    dplyr::n_distinct(reads$sample_id), " samples\n", sep = "")

reads <- reads |>
  dplyr::mutate(any_soft = (softL > 0) | (softR > 0),
                pos_rel  = pos - CDS_END)   # bp from 3' edge of CDS window

# ---- 2. Classify samples by fbxl1 introgression -------------------------
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

sample_df <- metadata |>
  dplyr::filter(sample_id %in% reads$sample_id, plate %in% c(1, 2, 3, 4)) |>
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
cat("\nfbxl1 introgression classification:\n"); print(table(sample_df$Genotype))

# ---- 3. Join and per-sample summary --------------------------------------
reads <- reads |> dplyr::inner_join(
  sample_df |> dplyr::select(sample_id, Genotype, taxa), by = "sample_id")

per_sample <- reads |>
  dplyr::group_by(sample_id, Genotype, taxa) |>
  dplyr::summarise(
    n_reads         = dplyr::n(),
    n_softclipped   = sum(any_soft),
    softclip_rate   = mean(any_soft),
    mean_softR      = mean(softR),
    frac_softR_gt5  = mean(softR > 5),
    .groups = "drop"
  )

cat("\nPer-group summary:\n")
per_sample |> dplyr::group_by(Genotype) |>
  dplyr::summarise(
    n_samples          = dplyr::n(),
    total_reads        = sum(n_reads),
    mean_n_reads       = mean(n_reads),
    mean_softclip_rate = mean(softclip_rate),
    mean_softR         = mean(mean_softR),
    .groups = "drop"
  ) |> print()

write.csv(per_sample, file.path(out_dir, "per_sample_summary.csv"),
          row.names = FALSE)

# ---- 4. Plots ------------------------------------------------------------
col_geno <- c("B73" = "#03bec4", "Teo" = "#d62728")

# soft-clip rate per sample, split by group
p_rate <- ggplot(per_sample |> dplyr::filter(n_reads >= 20),
                 aes(x = Genotype, y = softclip_rate, fill = Genotype)) +
  geom_boxplot(width = 0.4, outlier.size = 0.7, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.9) +
  scale_fill_manual(values = col_geno, guide = "none") +
  labs(
    title    = "Soft-clip rate in CDS window, per sample",
    subtitle = "Fraction of reads carrying any soft-clip. Samples with <20 reads dropped.",
    y = "fraction reads soft-clipped", x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "softclip_rate_by_group.png"),
       p_rate, width = 6, height = 5, dpi = 200)

# start-position distribution across the window
# samtools view returns reads that OVERLAP the region, so leftmost POS can be
# far upstream for spliced reads (large N in CIGAR). Keep only reads whose
# start actually falls within the window itself for the density.
reads_in_window <- reads |>
  dplyr::filter(pos >= CDS_START, pos <= CDS_END)

p_start <- ggplot(reads_in_window, aes(x = pos_rel, colour = Genotype)) +
  geom_density(linewidth = 0.7) +
  geom_vline(xintercept = c(CDS_START - CDS_END, 0),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = col_geno) +
  scale_x_continuous(limits = c(CDS_START - CDS_END, 10)) +
  labs(
    title    = "Read start-position distribution within CDS window",
    subtitle = sprintf("Reads whose 5' end falls inside the window (%d B73, %d Teo). Displaced UTR reads cluster near 0 (3' edge).",
                       sum(reads_in_window$Genotype == "B73"),
                       sum(reads_in_window$Genotype == "Teo")),
    x = "bp from 3' edge of CDS window (0 = 17933902)",
    y = "density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "start_position_density.png"),
       p_start, width = 8, height = 5, dpi = 200)

# trailing soft-clip length density (3'-edge signature)
p_softR <- ggplot(reads |> dplyr::filter(softR > 0),
                  aes(x = softR, colour = Genotype)) +
  geom_density(linewidth = 0.7) +
  scale_colour_manual(values = col_geno) +
  labs(
    title    = "Trailing soft-clip length distribution (softR > 0)",
    subtitle = "Reads whose 3' end was clipped. Higher, longer clips = more UTR spill-in.",
    x = "trailing soft-clip length (bp)",
    y = "density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "softR_length_density.png"),
       p_softR, width = 8, height = 5, dpi = 200)

# ---- 5. Per-taxon stratification (Rubén §9 — check the distal taxa) -----
# Teo carriers only; sample-level means so a chatty sample can't dominate.
teo_per_sample <- per_sample |> dplyr::filter(Genotype == "Teo")

per_taxon <- teo_per_sample |>
  dplyr::group_by(taxa) |>
  dplyr::summarise(
    n_samples          = dplyr::n(),
    total_reads        = sum(n_reads),
    mean_softclip_rate = mean(softclip_rate),
    mean_softR         = mean(mean_softR),
    frac_softR_gt5_mean = mean(frac_softR_gt5),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(mean_softclip_rate))

# B73 reference row (from all B73 samples)
b73_ref <- per_sample |> dplyr::filter(Genotype == "B73") |>
  dplyr::summarise(
    taxa               = "B73_ref",
    n_samples          = dplyr::n(),
    total_reads        = sum(n_reads),
    mean_softclip_rate = mean(softclip_rate),
    mean_softR         = mean(mean_softR),
    frac_softR_gt5_mean = mean(frac_softR_gt5)
  )
per_taxon_full <- dplyr::bind_rows(per_taxon, b73_ref)

cat("\nPer-taxon soft-clip summary (Teo carriers, plus B73 reference):\n")
print(per_taxon_full, n = Inf)
write.csv(per_taxon_full,
          file.path(out_dir, "per_taxon_summary.csv"),
          row.names = FALSE)

# Flag: taxa Rubén hasn't verified for CDS conservation
unverified_taxa <- c("Zlux", "Zdip", "Hueh")

# bar plot: soft-clip rate per taxon, with B73 reference line
b73_rate <- b73_ref$mean_softclip_rate
b73_softR <- b73_ref$mean_softR

p_tax_rate <- ggplot(per_taxon,
                     aes(x = reorder(taxa, mean_softclip_rate),
                         y = mean_softclip_rate,
                         fill = taxa %in% unverified_taxa)) +
  geom_hline(yintercept = b73_rate, linetype = "dashed", colour = "grey40") +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = paste0("n=", n_samples)),
            vjust = -0.3, size = 3.4) +
  scale_fill_manual(
    values = c(`TRUE` = "#f4a261", `FALSE` = "#03bec4"),
    labels = c(`TRUE` = "CDS conservation NOT verified (Rubén §9)",
               `FALSE` = "CDS conservation verified"),
    name   = NULL) +
  labs(
    title    = "Soft-clip rate per Teo donor taxon (CDS window)",
    subtitle = sprintf("Dashed line: B73 background rate = %.1f%%. Orange bars: distal taxa Rubén flagged as unverified.",
                       100 * b73_rate),
    x = NULL, y = "mean fraction reads soft-clipped"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top")

ggsave(file.path(out_dir, "softclip_rate_by_taxon.png"),
       p_tax_rate, width = 8, height = 5, dpi = 200)

# same for mean trailing clip length
p_tax_softR <- ggplot(per_taxon,
                      aes(x = reorder(taxa, mean_softR),
                          y = mean_softR,
                          fill = taxa %in% unverified_taxa)) +
  geom_hline(yintercept = b73_softR, linetype = "dashed", colour = "grey40") +
  geom_col(width = 0.65, colour = "black") +
  geom_text(aes(label = paste0("n=", n_samples)),
            vjust = -0.3, size = 3.4) +
  scale_fill_manual(
    values = c(`TRUE` = "#f4a261", `FALSE` = "#03bec4"),
    labels = c(`TRUE` = "CDS conservation NOT verified",
               `FALSE` = "CDS conservation verified"),
    name   = NULL) +
  labs(
    title    = "Mean trailing soft-clip length per Teo donor taxon",
    subtitle = sprintf("Dashed line: B73 background = %.2f bp. Bigger bars = more sequence divergence at the read's 3' end.",
                       b73_softR),
    x = NULL, y = "mean softR (bp)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top")

ggsave(file.path(out_dir, "softR_length_by_taxon.png"),
       p_tax_softR, width = 8, height = 5, dpi = 200)

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "per_sample_summary.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "per_taxon_summary.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "softclip_rate_by_group.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "start_position_density.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "softR_length_density.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "softclip_rate_by_taxon.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "softR_length_by_taxon.png"), "\n", sep = "")
