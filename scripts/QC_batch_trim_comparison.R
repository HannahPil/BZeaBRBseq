#!/usr/bin/env Rscript
# ==============================================================================
# QC -- Sequencing depth and Trimmomatic read-pair fate per plate, batch 1 vs 2
#
# Input:  data/trim_stats_by_plate.csv, transcribed from the "Input Read Pairs /
#         Both Surviving / ..." line of each trimmed/pool_N_trim.log on Hazel.
#         Add a row when another plate finishes trimming and rerun.
# Output: output/QC/batch_trim_comparison.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

d <- read.csv("data/trim_stats_by_plate.csv")
d$plate_lab <- factor(paste("Plate", d$plate), levels = paste("Plate", sort(d$plate)))
d$batch_lab <- factor(paste("Batch", d$batch))
batch_cols  <- c("Batch 1" = "#7F7F7F", "Batch 2" = "#0072B2")

# --- A: read pairs per plate ---------------------------------------------------
d$per_well <- d$input_pairs / 96
pA <- ggplot(d, aes(plate_lab, input_pairs / 1e6, fill = batch_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f M\nper well", per_well / 1e6)),
            vjust = -0.2, size = 4.6, lineheight = 0.9) +
  facet_grid(. ~ batch_lab, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = batch_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Read pairs (millions)", tag = "A")

# --- B: read-pair fate after trimming ------------------------------------------
fate_levels <- c("R1 and R2 kept", "R1 only kept", "R2 only kept", "R1 and R2 dropped")
long <- rbind(
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = "R1 and R2 kept", n = d$both_surviving),
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = "R1 only kept",    n = d$forward_only),
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = "R2 only kept",    n = d$reverse_only),
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = "R1 and R2 dropped",    n = d$dropped))
long$fate <- factor(long$fate, levels = rev(fate_levels))
long$pct  <- 100 * long$n / d$input_pairs[match(long$plate_lab, d$plate_lab)]

pB <- ggplot(long, aes(plate_lab, pct, fill = fate)) +
  geom_col(width = 0.7) +
  geom_text(data = long[long$fate == "R1 and R2 kept", ],
            aes(y = pct / 2, label = sprintf("%.1f%%", pct)),
            colour = "white", size = 4.6, fontface = "bold") +
  scale_fill_manual(values = c("R1 and R2 kept" = "#0072B2", "R1 only kept" = "#56B4E9",
                               "R2 only kept" = "#E69F00", "R1 and R2 dropped" = "#B3252B"),
                    breaks = fate_levels, name = NULL) +
  facet_grid(. ~ batch_lab, scales = "free_x", space = "free_x", switch = "x") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Read pairs (%)", tag = "B")

theme_set(theme_bw(base_size = 17) +
          theme(panel.background = element_rect(fill = "#ffffff"),
                plot.background  = element_rect(fill = "#ffffff", colour = NA),
                panel.grid.major.x = element_blank(),
                panel.grid.minor   = element_blank(),
                legend.position = "top",
                plot.title.position = "plot",
                plot.caption.position = "plot",
                plot.caption = element_text(hjust = 0, size = 13),
                plot.subtitle = element_text(size = 14),
                strip.placement = "outside",
                strip.background = element_blank(),
                strip.text = element_text(size = 17, face = "bold")))

missing <- setdiff(c(1:5, 7:9), d$plate)
caption <- paste0(
  sprintf("Batch 2 plates received %.1f-%.1f times as many read pairs as batch 1 plates.\n",
          min(d$input_pairs[d$batch == 2]) / max(d$input_pairs[d$batch == 1]),
          max(d$input_pairs[d$batch == 2]) / min(d$input_pairs[d$batch == 1])),
  "Pairs keeping R1 and R2 after trimming: batch 1 ",
  sprintf("%.1f-%.1f%%", min(100 * d$both_surviving[d$batch == 1] / d$input_pairs[d$batch == 1]),
                         max(100 * d$both_surviving[d$batch == 1] / d$input_pairs[d$batch == 1])),
  ", batch 2 ",
  sprintf("%.1f-%.1f%%", min(100 * d$both_surviving[d$batch == 2] / d$input_pairs[d$batch == 2]),
                         max(100 * d$both_surviving[d$batch == 2] / d$input_pairs[d$batch == 2])),
  ". Only pairs that keep R1 and R2 enter STARsolo.\n",
  "Plates 1-3 and 7 lost the largest shares of pairs to removal of R1 (the barcode and UMI read; 'R2 only kept').",
  if (length(missing)) paste0("\nPlate ", paste(missing, collapse = ", "),
                              ": trimming not yet complete.") else "")

p <- (pA / pB) + plot_layout(guides = "keep") +
  plot_annotation(
    title    = "BRB-seq read pairs and Trimmomatic read-pair fate per plate, batches 1 and 2",
    subtitle = paste0("Batch 1: NovaSeq 6000 S4, 1 lane, plates 1-4. Batch 2: NovaSeq X Plus 25B, 1 lane, plates 5, 7, 8, 9.\n",
                      "Trimmomatic PE (Nextera adapters, SLIDINGWINDOW 4:15, MINLEN 28); per well = read pairs / 96."),
    caption  = caption)

dir.create("output/QC", showWarnings = FALSE, recursive = TRUE)
ggsave("output/QC/batch_trim_comparison.png", p, width = 18, height = 12, dpi = 200, bg = "#ffffff")
cat("Wrote output/QC/batch_trim_comparison.png\n")
