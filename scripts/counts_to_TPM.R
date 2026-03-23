#!/usr/bin/env Rscript
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
data_dir    <- "data"
gtf_file    <- file.path(data_dir, "Zea_mays.gtf")
counts_file <- file.path(data_dir, "Zea_mays_counts.txt")
output_file <- file.path("output", "Zea_mays_TPM.txt")

# ── 1. Load counts ─────────────────────────────────────────────────────────────
counts <- read.table(counts_file, header = TRUE, sep = "\t", 
                     row.names = 1, check.names = FALSE)

cat("Loaded counts:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ── 2. Parse GTF to get gene lengths (non-overlapping exon bases) ──────────────
cat("Parsing GTF...\n")

gtf <- read.table(gtf_file, sep = "\t", quote = "", comment.char = "#",
                  col.names = c("chr", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attributes"))

# Keep only exon rows
exons <- gtf[gtf$feature == "exon", ]

# Extract gene_id
exons$gene_id <- gsub('.*gene_id "([^"]+)".*', "\\1", exons$attributes)

# Calculate non-overlapping length per gene by merging overlapping exons
gene_lengths <- exons %>%
  group_by(gene_id) %>%
  arrange(start, .by_group = TRUE) %>%
  summarise(length = {
    starts <- start
    ends   <- end
    merged_start <- starts[1]
    merged_end   <- ends[1]
    total_len    <- 0
    for (i in seq_along(starts)) {
      if (i == 1) next
      if (starts[i] <= merged_end) {
        merged_end <- max(merged_end, ends[i])
      } else {
        total_len  <- total_len + (merged_end - merged_start + 1)
        merged_start <- starts[i]
        merged_end   <- ends[i]
      }
    }
    total_len + (merged_end - merged_start + 1)
  }, .groups = "drop")

cat("Gene lengths calculated for", nrow(gene_lengths), "genes\n")

# ── 3. Align genes between counts and GTF ─────────────────────────────────────
common_genes <- intersect(rownames(counts), gene_lengths$gene_id)
cat("Genes in common:", length(common_genes), "\n")
cat("Genes in counts but not GTF:", 
    length(setdiff(rownames(counts), gene_lengths$gene_id)), "\n")

counts       <- counts[common_genes, ]
gene_lengths <- gene_lengths[match(common_genes, gene_lengths$gene_id), ]

# ── 4. Calculate TPM ───────────────────────────────────────────────────────────
len_kb <- gene_lengths$length / 1000
rpk    <- counts / len_kb
tpm    <- t(t(rpk) / (colSums(rpk) / 1e6))

cat("TPM calculation complete. Column sums (should all be ~1,000,000):\n")
print(round(colSums(tpm)))

# ── 5. Save output ─────────────────────────────────────────────────────────────
write.table(tpm, file = output_file, sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = NA)

cat("TPM saved to:", output_file, "\n")