# ============================
# setup
# ============================
install.packages("ggplot2")
install.packages("ggrepel")

library(edgeR)
library(ggplot2)
library(ggrepel)

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea RNA-seq/BZeaBRBseq")

# ============================
# load data
# ============================
counts <- read.table("Zea_mays_counts.txt", header = TRUE, row.names = 1, sep = "\t")
meta_all <- read.table("metadata.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")

# keep only complete rows (plate 1–4)
meta <- subset(meta_all, plate %in% c(1, 2, 3, 4))

# mark checks (B73 vs nonB73)
meta$check_status <- ifelse(meta$genotype == "B73", "check", "noncheck")

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
# plot
# ============================
p <- ggplot(scores, aes(x = PC1, y = PC2, color = check_status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_bw() +
  theme(legend.title = element_blank())

print(p)
