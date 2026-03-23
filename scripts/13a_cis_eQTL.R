# ==============================================================================
# CIS-EQTL ANALYSIS PIPELINE
# ==============================================================================
# Analysis of teosinte introgression effects on gene expression in maize
# Using MatrixEQTL for association testing
# Model: expression ~ introgression_status + plate
# Outputs: saved to ./output/cis_eQTL
# ==============================================================================

library(tidyverse)
library(MatrixEQTL)
library(rtracklayer)
library(GenomicRanges)

# ==============================================================================
# 0. OUTPUT DIR
# ==============================================================================

out_dir <- file.path("output", "cis_eQTL")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

print("1. Loading data...")

data_dir <- "data"
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
counts   <- read.delim(file.path(data_dir, "Zea_mays_counts.txt"), check.names = FALSE, row.names = 1)
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)

teogeno <- teogeno[!duplicated(names(teogeno))]

# ==============================================================================
# 2. ALIGN SAMPLES
# ==============================================================================

print("2. Aligning samples...")

sample_df <- metadata |>
  filter(
    sample_id %in% colnames(counts),
    plate %in% c(1, 2, 3, 4)                # keep only sequenced plates
  ) |>
  transmute(
    sample_id,
    taxa  = factor(taxa),
    plate = factor(plate),
    genotype_teogeno_key = if_else(genotype == "B73", "B73.B", paste0(genotype, ".B"))
  ) |>
  filter(genotype_teogeno_key %in% names(teogeno)) |>
  arrange(sample_id)

counts <- counts[, as.character(sample_df$sample_id)]

# ==============================================================================
# 3. PROCESS GTF COORDINATES
# ==============================================================================

print("3. Processing GTF coordinates...")

gtf <- import(file.path(data_dir, "Zea_mays.gtf"))
gtf_df <- as.data.frame(gtf)

gene_coords <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(
    chr = as.character(seqnames[1]),
    start = min(start),
    end   = max(end),
    .groups = "drop"
  ) |>
  filter(gene_id %in% rownames(counts)) |>
  arrange(gene_id)

genes_gr <- GRanges(
  seqnames = gene_coords$chr,
  ranges   = IRanges(gene_coords$start, gene_coords$end),
  gene_id  = gene_coords$gene_id
)

# ==============================================================================
# 4. BUILD GENOTYPE MATRIX FROM INTROGRESSIONS
# ==============================================================================

print("4. Building genotype matrix from introgressions...")

seg_df <- imap_dfr(teogeno, ~{
  .x |>
    filter(V4 == "Introgression") |>
    transmute(
      key = .y,
      chr = V1,
      start = as.integer(V2),
      end   = as.integer(V3)
    )
})

segs_gr <- GRanges(
  seqnames = seg_df$chr,
  ranges   = IRanges(seg_df$start, seg_df$end),
  key      = seg_df$key
)

hits <- findOverlaps(segs_gr, genes_gr)

genotype_mat <- matrix(
  0,
  nrow = nrow(gene_coords),
  ncol = nrow(sample_df),
  dimnames = list(gene_coords$gene_id, sample_df$sample_id)
)

hit_genes <- genes_gr$gene_id[subjectHits(hits)]
hit_keys  <- segs_gr$key[queryHits(hits)]
key_to_samples <- split(sample_df$sample_id, sample_df$genotype_teogeno_key)

for(i in seq_along(hit_genes)) {
  g_id <- hit_genes[i]
  t_key <- hit_keys[i]
  s_ids <- key_to_samples[[t_key]]
  if(!is.null(s_ids)) {
    genotype_mat[g_id, as.character(s_ids)] <- 1
  }
}

# ==============================================================================
# 5. PREPARE MATRIXEQTL INPUTS
# ==============================================================================

print("5. Preparing MatrixEQTL inputs...")

var_mask <- rowSums(genotype_mat) > 0 & rowSums(genotype_mat) < ncol(genotype_mat)
genotype_mat <- genotype_mat[var_mask, , drop = FALSE]

final_genes <- rownames(genotype_mat)
gene_coords_subset <- gene_coords[match(final_genes, gene_coords$gene_id), ]

# NOTE: using raw library-size CPM (not TMM-normalized) here intentionally.
# MatrixEQTL's linear model with plate covariates handles composition differences.
# The PCA/DE scripts (10, 11, 12) use edgeR TMM CPM instead.
lib_size <- colSums(counts)
expr_mat <- as.matrix(counts[rownames(genotype_mat), ])
expr_mat <- log2((t(t(expr_mat) / lib_size)) * 1e6 + 1)

gene_slice <- SlicedData$new()
gene_slice$CreateFromMatrix(expr_mat)
gene_slice$ResliceCombined(sliceSize = 2000)

snps_slice <- SlicedData$new()
snps_slice$CreateFromMatrix(genotype_mat)
snps_slice$ResliceCombined(sliceSize = 2000)

print("   > covariates: plate only")
cvrt_mat <- model.matrix(~ plate, data = sample_df)
cvrt_mat <- cvrt_mat[, -1, drop = FALSE]
cvrt_mat <- t(cvrt_mat)

print(paste("   > covariates used:", paste(rownames(cvrt_mat), collapse = ", ")))

cvrt_slice <- SlicedData$new()
cvrt_slice$CreateFromMatrix(cvrt_mat)

snpspos <- as.data.frame(gene_coords_subset[, c("gene_id", "chr", "start")])
colnames(snpspos) <- c("snp", "chr", "pos")

genepos <- as.data.frame(gene_coords_subset[, c("gene_id", "chr", "start", "end")])
colnames(genepos) <- c("geneid", "chr", "left", "right")

# ==============================================================================
# 6. RUN MATRIXEQTL
# ==============================================================================

print("6. Running MatrixEQTL...")

me <- Matrix_eQTL_main(
  snps = snps_slice,
  gene = gene_slice,
  cvrt = cvrt_slice,
  output_file_name = NULL,
  pvOutputThreshold = 0,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = file.path(out_dir, "BZea_cis_eQTL_results.csv"),
  pvOutputThreshold.cis = 1,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 10,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

print(paste("SUCCESS! Results saved to", file.path(out_dir, "BZea_cis_eQTL_results.csv")))

# ==============================================================================
# 7. LOAD AND PREPARE RESULTS
# ==============================================================================

print("7. Loading results...")

results <- read.delim(
  file.path(out_dir, "BZea_cis_eQTL_results.csv"),
  sep = "\t",
  stringsAsFactors = FALSE
)

write.csv(results, file.path(out_dir, "BZea_cis_eQTL_results_excel.csv"), row.names = FALSE)

n_sig <- sum(results$FDR < 0.05)
print(paste("Significant hits (FDR < 0.05):", n_sig))
print(paste("Total tests:", nrow(results)))

results <- results |>
  mutate(
    is_sig = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    log_p = -log10(p.value)
  )

sig_results <- results |> filter(FDR < 0.05)
write.csv(sig_results, file.path(out_dir, "BZea_Significant_eQTLs.csv"), row.names = FALSE)

top_hits <- results |>
  arrange(p.value) |>
  head(10)

print("Top 10 hits:")
print(top_hits)

# ==============================================================================
# 8. GENERATE PLOTS
# ==============================================================================

print("8. Generating plots...")

p1 <- ggplot(results, aes(x = p.value)) +
  geom_histogram(bins = 50, fill = "dodgerblue", color = "white") +
  theme_minimal() +
  labs(
    title = "P-value distribution",
    subtitle = paste("Total tests:", nrow(results)),
    x = "p-value",
    y = "frequency"
  )

ggsave(file.path(out_dir, "plot_pvalue_hist.png"), p1, width = 6, height = 4)

p2 <- ggplot(results, aes(x = beta, y = log_p, color = is_sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Volcano plot: introgression effects",
    subtitle = paste(n_sig, "significant genes (FDR < 0.05)"),
    x = "effect size (beta)",
    y = "-log10(p-value)",
    color = "status"
  ) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "plot_volcano.png"), p2, width = 6, height = 5)

plot_data <- results |>
  inner_join(snpspos, by = c("SNP" = "snp")) |>
  mutate(
    chr_clean = gsub("chr", "", chr, ignore.case = TRUE),
    chr_num = as.numeric(chr_clean),
    log_p = -log10(p.value)
  ) |>
  filter(!is.na(chr_num))

data_cum <- plot_data |>
  group_by(chr_num) |>
  summarise(max_bp = as.numeric(max(pos)), .groups = "drop") |>
  arrange(chr_num) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr_num, bp_add)

plot_data <- plot_data |>
  inner_join(data_cum, by = "chr_num") |>
  mutate(
    bp_cum = as.numeric(pos + bp_add),
    color_group = case_when(
      FDR < 0.05 & (chr_num %% 2 == 1) ~ "Sig_Odd",
      FDR < 0.05 & (chr_num %% 2 == 0) ~ "Sig_Even",
      (chr_num %% 2 == 1) ~ "Base_Odd",
      TRUE ~ "Base_Even"
    )
  ) |>
  arrange(FDR)

axis_set <- plot_data |>
  group_by(chr_num) |>
  summarize(
    center = (as.numeric(max(bp_cum)) + as.numeric(min(bp_cum))) / 2,
    .groups = "drop"
  ) |>
  arrange(chr_num)

sig_cutoff <- -log10(max(results$p.value[results$FDR < 0.05]))
max_y <- max(plot_data$log_p, na.rm = TRUE)

p_man <- ggplot(plot_data, aes(x = bp_cum, y = log_p, color = color_group)) +
  geom_point(alpha = 0.75, size = 1.3) +
  scale_color_manual(values = c(
    "Base_Odd"  = "grey45",
    "Base_Even" = "grey70",
    "Sig_Odd"   = "firebrick",
    "Sig_Even"  = "darkred"
  )) +
  geom_hline(yintercept = sig_cutoff, color = "red", linetype = "dashed") +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = as.character(axis_set$chr_num)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max_y * 1.1)
  ) +
  labs(
    title = "Cis-eQTL manhattan plot",
    subtitle = paste("Genome-wide view of", nrow(plot_data), "genes"),
    x = "chromosome",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 10, color = "black")
  )

ggsave(file.path(out_dir, "ciseqtl_manhattan.png"), p_man, width = 10, height = 5)

