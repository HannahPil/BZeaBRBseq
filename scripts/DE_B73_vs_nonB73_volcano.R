suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
})

cts <- as.matrix(read.table("data/processed/Zea_mays_counts.txt",
                            header = TRUE, row.names = 1, check.names = FALSE,
                            sep = "\t"))
meta <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)
meta <- meta[match(colnames(cts), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(cts)))

# Keep only B73 checks (check == "yes") and introgression lines. Drop wells
# that are neither, which would otherwise fall into "nonB73":
#   - blanks (check == "blank", added with batch 2)
#   - Purple Checks (check == "no" in metadata, but not introgression lines)
#   - wells with no genotype recorded
is_line <- meta$check == "no" & !is.na(meta$genotype) & meta$genotype != "NA" &
           trimws(meta$genotype) != "Purple Check"
drop <- !(meta$check == "yes" | is_line)
if (any(drop)) {
  cat("Dropping", sum(drop), "non-B73, non-line well(s):\n")
  print(table(ifelse(meta$check[drop] == "blank", "blank", meta$genotype[drop]), useNA = "ifany"))
  cts  <- cts[, !drop, drop = FALSE]
  meta <- meta[!drop, ]
}

meta$group <- factor(ifelse(meta$check == "yes", "B73", "nonB73"),
                     levels = c("B73", "nonB73"))
meta$plate <- factor(meta$plate)

keep <- rowSums(cts >= 5) >= 10
cts  <- cts[keep, ]

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = meta,
                              design    = ~ plate + group)
dds <- DESeq(dds, parallel = FALSE)
res <- results(dds, name = "group_nonB73_vs_B73",
               independentFiltering = TRUE)

out_tab <- as.data.frame(res)
out_tab$gene <- rownames(out_tab)
out_tab <- out_tab[order(out_tab$padj), ]
dir.create("output/B73_vs_nonB73", showWarnings = FALSE, recursive = TRUE)
fwrite(out_tab, "output/B73_vs_nonB73/DE_results.csv")

df <- out_tab[!is.na(out_tab$padj), ]
ymax <- 13
df$neglog10p_clip <- pmin(-log10(df$padj), ymax)
df$sig <- df$padj < 0.05 & abs(df$log2FoldChange) >= 1
df$plotgrp <- ifelse(df$sig,
                     ifelse(df$log2FoldChange > 0,
                            "Higher in non-B73 (FDR<0.05, |log2FC|≥1)",
                            "Higher in B73 (FDR<0.05, |log2FC|≥1)"),
                     "ns")

label_ppl2 <- df[df$gene == "Zm00001eb012750", ]
xlim <- max(abs(df$log2FoldChange), na.rm = TRUE)

p <- ggplot(df, aes(log2FoldChange, neglog10p_clip, color = plotgrp)) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  scale_color_manual(values = c(
    "Higher in non-B73 (FDR<0.05, |log2FC|≥1)" = "#D7261E",
    "Higher in B73 (FDR<0.05, |log2FC|≥1)"     = "#1F6FEB",
    "ns" = "grey75"), name = NULL) +
  { if (nrow(label_ppl2) == 1)
      geom_text_repel(data = label_ppl2,
        aes(label = "ppl2 (Zm00001eb012750)"),
        color = "black", size = 3.3, min.segment.length = 0,
        nudge_x = 1.5, nudge_y = 1, box.padding = 0.5,
        show.legend = FALSE, inherit.aes = TRUE)
  } +
  scale_x_continuous(limits = c(-xlim, xlim)) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Differential expression: non-B73 (n = 348) vs B73 (n = 36)",
    subtitle = sprintf(
      "DESeq2 Wald, design ~ plate + group; %d genes tested; dashed line = FDR 0.05, |log2FC|≥1",
      nrow(df)),
    x = "log2 fold change (B73×teosinte / B73)",
    y = "-log10 adjusted p",
    caption = "348 teosinte-donor BRB-seq lines vs 36 B73 checks on the same 4 plates. Genes to the right are up in lines carrying teosinte DNA."
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        legend.direction = "vertical",
        plot.title.position = "plot",
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.9, color = "black"),
        axis.ticks = element_line(linewidth = 0.7, color = "black"))

ggsave("output/B73_vs_nonB73/volcano_B73_vs_nonB73.png",
       p, width = 5, height = 6.5, dpi = 200)
