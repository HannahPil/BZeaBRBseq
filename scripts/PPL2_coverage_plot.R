#!/usr/bin/env Rscript

# ==============================================================================
# PPL2 coverage (LOCAL) — per-base coverage profile across ppl2 / PnsL1
#
# Decides between the two explanations left standing after the in-silico
# artifact checks (BZea_PPL2_2026/scripts/11):
#
#   (a) B73 genuinely does not express ppl2
#          -> B73 trace flat across the whole locus
#   (b) B73 expresses it but the 3' gene model is wrong, so BRB-seq counts
#       the wrong window
#          -> B73 trace has signal over the gene body, absent only at the 3' end
#
# A third pattern worth watching for: signal in B73 that sits PAST the
# annotated 3' end. ppl2 (+) and Zm00001eb012760 (-) are convergent and only
# 25 bp apart, so anything downstream is more likely the neighbour's 3' end
# than an extended ppl2 UTR -- which is why the window includes both flanks.
#
# INPUT (produced on HPC by PPL2_coverage_hpc.sh, shipped via git):
#   data/PPL2_depth_matrix.tsv
#
# OUTPUT:
#   output/PPL2_coverage/ppl2_coverage_profile.png
#   output/PPL2_coverage/ppl2_coverage_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

FOCUS_GENE   <- "Zm00001eb012750"   # ppl2 / PnsL1
GENE_START   <- 42030859L; GENE_END <- 42032683L    # + strand, 1825 bp
NBR_UP_END   <- 42030276L           # Zm00001eb012740 ends here
NBR_DN       <- "Zm00001eb012760"
NBR_DN_START <- 42032708L           # - strand, so this end is its 3' terminus
NBR_DN_END   <- 42041803L           # ... and this end is its 5'; 9.1 kb total
# The ppl2 panel keeps the original view; the depth matrix now runs to 42042300
# so the neighbour can be profiled over its whole length in section 6.
VIEW_START   <- 42029800L; VIEW_END <- 42033600L

data_dir <- "data"
out_dir  <- file.path("output", "PPL2_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Depth matrix ------------------------------------------------------
depth <- read.table(file.path(data_dir, "PPL2_depth_matrix.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE)
stopifnot(all(c("chr", "pos") %in% colnames(depth)))
cat("Depth matrix: ", nrow(depth), " positions x ", ncol(depth) - 2,
    " samples\n", sep = "")

# ---- 2. Library sizes, from the same matrix the checks used ---------------
# Normalising by library size so the two groups are comparable; the counts
# matrix is the canonical UMI-collapsed one from the Sep 2026 rerun.
counts <- read.delim(file.path(data_dir, "processed", "Zea_mays_counts.txt"),
                     check.names = FALSE, row.names = 1)
lib_vec <- colSums(counts)

# ---- 3. Classify each sample by introgression AT ppl2 ---------------------
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

sample_df <- metadata |>
  dplyr::filter(plate %in% c(1, 2, 3, 4)) |>
  dplyr::transmute(sample_id, genotype,
    key = dplyr::if_else(genotype == "B73", "B73.B", paste0(genotype, ".B"))) |>
  dplyr::filter(key %in% names(teogeno))

seg_df <- purrr::imap_dfr(teogeno, ~{
  .x |> dplyr::filter(V4 == "Introgression") |>
    dplyr::transmute(key = .y, chr = V1, start = as.integer(V2),
                     end = as.integer(V3))
})
segs_gr  <- GRanges(seg_df$chr, IRanges(seg_df$start, seg_df$end), key = seg_df$key)
focus_gr <- GRanges("chr1", IRanges(GENE_START, GENE_END))
carrier_keys <- unique(segs_gr$key[queryHits(findOverlaps(segs_gr, focus_gr))])

sample_df <- sample_df |>
  dplyr::mutate(group = dplyr::if_else(key %in% carrier_keys,
                                       "Teosinte carrier", "B73 background"))

samples <- intersect(sample_df$sample_id,
                     intersect(colnames(depth), names(lib_vec)))
sample_df <- sample_df |> dplyr::filter(sample_id %in% samples)
cat("Usable samples: ", nrow(sample_df), " (",
    sum(sample_df$group == "Teosinte carrier"), " carriers, ",
    sum(sample_df$group == "B73 background"), " B73 background)\n", sep = "")

# ---- 4. Library-normalised mean depth per position, per group -------------
long <- depth |>
  dplyr::select(pos, dplyr::all_of(samples)) |>
  tidyr::pivot_longer(-pos, names_to = "sample_id", values_to = "depth") |>
  dplyr::left_join(sample_df |> dplyr::select(sample_id, group), by = "sample_id") |>
  dplyr::mutate(cpm = depth / lib_vec[sample_id] * 1e6)

prof <- long |>
  dplyr::group_by(group, pos) |>
  dplyr::summarise(mean_cpm = mean(cpm),
                   se_cpm   = sd(cpm) / sqrt(dplyr::n()),
                   .groups = "drop")

# exon blocks for the gene track underneath the profile
gtf <- import(file.path(data_dir, "external", "Zea_mays.gtf"))
exons <- as.data.frame(gtf) |>
  dplyr::filter(type == "exon", grepl(FOCUS_GENE, gene_id)) |>
  dplyr::distinct(start, end)

p <- ggplot(prof |> dplyr::filter(pos >= VIEW_START, pos <= VIEW_END),
            aes(pos, mean_cpm, colour = group, fill = group)) +
  annotate("rect", xmin = GENE_START, xmax = GENE_END,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.45) +
  annotate("rect", xmin = exons$start, xmax = exons$end,
           ymin = -Inf, ymax = 0, fill = "grey35", colour = NA) +
  geom_vline(xintercept = NBR_DN_START, linetype = "dashed",
             colour = "grey35", linewidth = 0.4) +
  annotate("text", x = NBR_DN_START, y = Inf, hjust = -0.05, vjust = 1.6,
           label = "Zm00001eb012760 (-) starts", size = 3.2, colour = "grey35") +
  geom_ribbon(aes(ymin = mean_cpm - se_cpm, ymax = mean_cpm + se_cpm),
              colour = NA, alpha = 0.22) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("B73 background" = "#1f77b4",
                                 "Teosinte carrier" = "#d62728")) +
  scale_fill_manual(values = c("B73 background" = "#1f77b4",
                               "Teosinte carrier" = "#d62728")) +
  labs(title = "Per-base coverage across ppl2 (Zm00001eb012750)",
       subtitle = paste("Grey block = annotated gene span; dark blocks = exons.",
                        "Mean +/- SE, library-normalised."),
       x = "chr1 position (bp)", y = "Depth (CPM)",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"),
        legend.position = "top")

print(p)
ggsave(file.path(out_dir, "ppl2_coverage_profile.png"), p,
       width = 10, height = 5.5, dpi = 300)

# ---- 5. Where does the signal sit? ---------------------------------------
# Splits the gene into a 3' window (what BRB-seq counts) and the rest. If (b)
# is true, B73 has body signal without 3' signal.
THREE_PRIME_FROM <- GENE_END - 300L
summ <- long |>
  dplyr::mutate(window = dplyr::case_when(
      pos <  GENE_START                        ~ "upstream of gene",
      pos >= GENE_START & pos < THREE_PRIME_FROM ~ "gene body (5' of last 300bp)",
      pos >= THREE_PRIME_FROM & pos <= GENE_END  ~ "3' 300 bp (BRB-seq counts here)",
      pos >  GENE_END                          ~ "downstream of gene")) |>
  dplyr::group_by(group, window) |>
  dplyr::summarise(mean_cpm = mean(cpm),
                   pct_samples_with_any = 100 * mean(depth > 0),
                   .groups = "drop")
cat("\n=== where the reads sit ===\n")
print(as.data.frame(summ), digits = 3, row.names = FALSE)
write.csv(summ, file.path(out_dir, "ppl2_coverage_summary.csv"), row.names = FALSE)

# ---- 6. the downstream neighbour, over its whole length -------------------
# Why: the first window reached only 893 bp into Zm00001eb012760 and carriers
# looked 2.16x HIGHER there, while gene-level counts put them at 0.71x. Either
# the stub was unrepresentative, or the excess at the stub is ppl2 read-through
# rather than neighbour signal. Profiling the full 9.1 kb separates them:
#   - deficit spread evenly across the body -> mapping loss on divergent
#     teosinte sequence, since the introgression covers this gene too
#   - deficit concentrated somewhere        -> regulatory
if (max(depth$pos) < NBR_DN_END) {
  cat("\n[section 6 skipped: depth matrix stops at ", max(depth$pos),
      ", before the neighbour's far end at ", NBR_DN_END,
      ". Rerun PPL2_coverage_hpc.sh with the widened REGION.]\n", sep = "")
} else {
  nbr <- long |>
    dplyr::filter(pos >= NBR_DN_START, pos <= NBR_DN_END) |>
    dplyr::mutate(offset = pos - GENE_END,          # distance past ppl2's 3' end
                  bin = (offset %/% 250L) * 250L) |>
    dplyr::group_by(bin, group) |>
    dplyr::summarise(mean_cpm = mean(cpm), .groups = "drop") |>
    tidyr::pivot_wider(names_from = group, values_from = mean_cpm) |>
    dplyr::rename(B73bg = `B73 background`, carrier = `Teosinte carrier`) |>
    dplyr::mutate(ratio = carrier / pmax(B73bg, 1e-9))

  cat("\n=== Zm00001eb012760 across its full 9.1 kb, 250 bp bins ===\n")
  cat("    offset = bp past ppl2's 3' end. Neighbour runs +25 to +9120.\n")
  print(as.data.frame(nbr |>
    dplyr::transmute(offset_bp = bin, B73bg = round(B73bg, 2),
                     carrier = round(carrier, 2), ratio = round(ratio, 2))),
    row.names = FALSE)
  write.csv(nbr, file.path(out_dir, "ppl2_neighbour_profile.csv"), row.names = FALSE)

  pn <- ggplot(nbr |> tidyr::pivot_longer(c(B73bg, carrier),
                                          names_to = "group", values_to = "cpm"),
               aes(bin, cpm, colour = group)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_colour_manual(values = c(B73bg = "#1f77b4", carrier = "#d62728"),
                        labels = c("B73 background", "Teosinte carrier")) +
    labs(title = paste0("Coverage across ", NBR_DN, " (the downstream neighbour)"),
         subtitle = paste("Minus strand: its 3' end is on the LEFT.",
                          "250 bp bins, library-normalised."),
         x = "bp past ppl2's 3' end", y = "Mean depth (CPM)", colour = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "top")
  print(pn)
  ggsave(file.path(out_dir, "ppl2_neighbour_profile.png"), pn,
         width = 10, height = 4.5, dpi = 300)
}

cat("\nSaved to ", out_dir, "\n", sep = "")
