

library(edgeR)
library(ggplot2)
library(ggrepel)

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea RNA-seq/BZeaBRBseq")

# ============================
# load data
# ============================
counts <- read.table("Zea_mays_counts.txt", header = TRUE, row.names = 1, sep = "\t")
meta_all <- read.csv("metadata.csv")

# keep only complete rows (plate 1–4)
meta <- subset(meta_all, plate %in% c(1, 2, 3, 4))

# restrict metadata to samples present in counts
meta <- meta[meta$sample_id %in% colnames(counts), ]

# reorder counts to match metadata order
counts <- counts[, meta$sample_id]

# ============================
# filter + normalize
# ============================
keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, ]

dge <- DGEList(counts = counts_filt)
dge <- calcNormFactors(dge)
mat <- cpm(dge, log = TRUE, prior.count = 1)

# ============================
# pca
# ============================
pca <- prcomp(t(mat))

scores <- as.data.frame(pca$x)
scores$sample_id <- rownames(scores)

# drop duplicate sample_id columns
scores <- scores[, !duplicated(names(scores))]

# merge metadata
scores <- merge(scores, meta, by = "sample_id")

# variance explained
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# ============================
# plotting (before correction)
# ============================

# by taxa (checks in black)
p1 <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = ifelse(check == "yes", "check", taxa)), size = 3) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "check" = "black",
      "Dura"  = "#1f77b4",
      "Nobo"  = "#ff7f0e",
      "Mesa"  = "#2ca02c",
      "Chal"  = "#d62728",
      "Bals"  = "#9467bd",
      "Zdip"  = "brown4",
      "Hueh"  = "yellow3",
      "Zlux"  = "#17becf"
    ),
    breaks = c("check", "Dura", "Nobo", "Mesa", "Chal", "Bals", "Zdip", "Hueh", "Zlux")
  ) +
  labs(color = "Taxa / Check",
       title = "PCA before batch correction (by taxa)")
p1

# by plate
p2 <- ggplot(scores, aes(x = PC1, y = PC2, color = as.factor(plate))) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_bw() +
  labs(color = "Plate",
       title = "PCA before batch correction (by plate)")
p2


# ============================
# plotting (after correction)
# ============================

# by taxa (checks in black)
p1_corr <- ggplot(scores_corr, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = ifelse(check == "yes", "check", taxa)), size = 3) +
  xlab(paste0("PC1 (", pct_corr[1], "%)")) +
  ylab(paste0("PC2 (", pct_corr[2], "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "check" = "black",
      "Dura"  = "#1f77b4",
      "Nobo"  = "#ff7f0e",
      "Mesa"  = "#2ca02c",
      "Chal"  = "#d62728",
      "Bals"  = "#9467bd",
      "Zdip"  = "brown4",
      "Hueh"  = "yellow3",
      "Zlux"  = "#17becf"
    ),
    breaks = c("check", "Dura", "Nobo", "Mesa", "Chal", "Bals", "Zdip", "Hueh", "Zlux")
  ) +
  labs(color = "Taxa / Check",
       title = "PCA after batch correction (by taxa)")
p1_corr

# by plate
p2_corr <- ggplot(scores_corr, aes(x = PC1, y = PC2, color = as.factor(plate))) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", pct_corr[1], "%)")) +
  ylab(paste0("PC2 (", pct_corr[2], "%)")) +
  theme_bw() +
  labs(color = "Plate",
       title = "PCA after batch correction (by plate)")
p2_corr

#export
ggsave("output/PCA_taxa.png", p1, width = 7, height = 7, dpi = 300)
ggsave("output/PCA_plate.png", p2, width = 7, height = 7, dpi = 300)
ggsave("output/PCA_corr_taxa.png", p1_corr, width = 7, height = 7, dpi = 300)
ggsave("output/PCA_corr_plate.png", p2_corr, width = 7, height = 7, dpi = 300)

# ========================================================================================
# subset to FT_genes
# ========================================================================================
FT_genes <- c("Zm00001eb391060", "Zm00001eb353250")

counts <- counts[rownames(counts) %in% FT_genes, ]

# ============================
# normalize
# ============================
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
mat <- cpm(dge, log = TRUE, prior.count = 1)

# ============================
# PCA
# ============================
pca <- prcomp(t(mat))
scores <- as.data.frame(pca$x)
scores$sample_id <- rownames(scores)

# merge metadata
scores <- merge(scores, meta, by = "sample_id")

# variance explained
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# ============================
# PCA plot (FT_genes only)
# ============================
p_FT <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = ifelse(check == "yes", "check", taxa)), size = 3) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "check" = "black",
      "Dura"  = "#1f77b4",
      "Nobo"  = "#ff7f0e",
      "Mesa"  = "#2ca02c",
      "Chal"  = "#d62728",
      "Bals"  = "#9467bd",
      "Zdip"  = "brown4",
      "Hueh"  = "yellow3",
      "Zlux"  = "#17becf"
    ),
    breaks = c("check", "Dura", "Nobo", "Mesa", "Chal", "Bals", "Zdip", "Hueh", "Zlux")
  ) +
  labs(color = "Taxa",
       title = "PCA using FT_genes only")
p_FT


# ============================
# export
# ============================
ggsave("output/PCA_FT_genes.png", p_FT, width = 7, height = 7, dpi = 300)

apply(mat, 1, var)

