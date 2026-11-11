#!/usr/bin/env Rscript
# ==============================================================================
# QC -- STARsolo mapping outcome per plate, batch 1 vs 2
#
# Input:  data/starsolo_mapping_by_plate.csv, transcribed from each
#         starsolo/pool_N/Log.final.out on Hazel ("Number of input reads",
#         "Uniquely mapped reads %", "% of reads mapped to too many loci").
#         input_reads = read pairs keeping R1 and R2 after trimming.
#         Add a row when another plate finishes STARsolo and rerun.
# Output: output/QC/batch_mapping_comparison.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

d <- read.csv("data/starsolo_mapping_by_plate.csv")
d <- d[order(d$plate), ]
d$plate_lab <- factor(paste("Plate", d$plate), levels = paste("Plate", d$plate))
d$batch_lab <- factor(paste("Batch", d$batch))
d$unmapped_pct <- 100 - d$unique_pct - d$too_many_loci_pct
d$unique_reads <- d$input_reads * d$unique_pct / 100
batch_cols <- c("Batch 1" = "#7F7F7F", "Batch 2" = "#0072B2")

# --- A: uniquely mapped reads per plate ----------------------------------------
pA <- ggplot(d, aes(plate_lab, unique_reads / 1e6, fill = batch_lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f M\nper well", unique_reads / 96 / 1e6)),
            vjust = -0.2, size = 4.6, lineheight = 0.9) +
  facet_grid(. ~ batch_lab, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = batch_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Uniquely mapped reads (millions)", tag = "A")

# --- B: mapping outcome -------------------------------------------------------
fate_levels <- c("Uniquely mapped", "Mapped to too many loci", "Unmapped (too short or other)")
long <- rbind(
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = fate_levels[1], pct = d$unique_pct),
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = fate_levels[2], pct = d$too_many_loci_pct),
  data.frame(plate_lab = d$plate_lab, batch_lab = d$batch_lab, fate = fate_levels[3], pct = d$unmapped_pct))
long$fate <- factor(long$fate, levels = rev(fate_levels))

lab <- long[long$fate != fate_levels[3], ]
lab$ypos <- ifelse(lab$fate == fate_levels[1], lab$pct / 2,
                   d$unique_pct[match(lab$plate_lab, d$plate_lab)] + lab$pct / 2)

pB <- ggplot(long, aes(plate_lab, pct, fill = fate)) +
  geom_col(width = 0.7) +
  geom_text(data = lab, aes(y = ypos, label = sprintf("%.1f%%", pct)),
            colour = "white", size = 4.6, fontface = "bold") +
  scale_fill_manual(values = c("Uniquely mapped" = "#0072B2", "Mapped to too many loci" = "#E69F00",
                               "Unmapped (too short or other)" = "#B3252B"),
                    breaks = fate_levels, name = NULL) +
  facet_grid(. ~ batch_lab, scales = "free_x", space = "free_x", switch = "x") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Reads (%)", tag = "B")

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

rng <- function(x) sprintf("%.1f-%.1f%%", min(x), max(x))
missing <- setdiff(c(1:5, 7:9), d$plate)
caption <- paste0(
  "Uniquely mapped: batch 1 ", rng(d$unique_pct[d$batch == 1]),
  ", batch 2 ", rng(d$unique_pct[d$batch == 2]), ". ",
  "Mapped to too many loci: batch 1 ", rng(d$too_many_loci_pct[d$batch == 1]),
  ", batch 2 ", rng(d$too_many_loci_pct[d$batch == 2]), ".\n",
  "Only uniquely mapped reads are counted; reads mapping to more than one locus are discarded.",
  if (length(missing)) paste0("\nPlate ", paste(missing, collapse = ", "),
                              ": STARsolo not yet complete.") else "")

p <- (pA / pB) +
  plot_annotation(
    title    = "STARsolo mapping outcome per plate, batches 1 and 2",
    subtitle = paste0("Input: read pairs keeping R1 and R2 after trimming. STARsolo, ",
                      "--outFilterMultimapNmax 1; per well = uniquely mapped reads / 96."),
    caption  = caption)

dir.create("output/QC", showWarnings = FALSE, recursive = TRUE)
ggsave("output/QC/batch_mapping_comparison.png", p, width = 18, height = 12, dpi = 200, bg = "#ffffff")
cat("Wrote output/QC/batch_mapping_comparison.png\n")
