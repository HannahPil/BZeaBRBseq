# ==============================================================================
# CIS-EQTL ANALYSIS PIPELINE
# ==============================================================================
# Analysis of teosinte introgression effects on gene expression in maize
# Using MatrixEQTL for association testing
# ==============================================================================

library(tidyverse)
library(MatrixEQTL)
library(rtracklayer)
library(GenomicRanges)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

print("1. Loading Data...")

# Load inputs
teogeno  <- readRDS("results_list_new_name.rds")
counts   <- read.delim("Zea_mays_counts.txt", check.names = FALSE, row.names = 1)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Clean duplicates from genotype list
teogeno <- teogeno[!duplicated(names(teogeno))]

# ==============================================================================
# 2. ALIGN SAMPLES
# ==============================================================================

print("2. Aligning Samples...")

sample_df <- metadata |>
  filter(sample_id %in% colnames(counts)) |>
  transmute(
    sample_id,
    taxa  = factor(taxa),
    plate = factor(plate),
    # Map Genotype to Teogeno Key (B73 -> B73.B)
    genotype_teogeno_key = if_else(genotype == "B73", "B73.B", paste0(genotype, ".B"))
  ) |>
  filter(genotype_teogeno_key %in% names(teogeno)) |>
  arrange(sample_id)

# Align counts to final sample list
counts <- counts[, as.character(sample_df$sample_id)]

# ==============================================================================
# 3. PROCESS GENE COORDINATES
# ==============================================================================

print("3. Processing GTF Coordinates...")

gtf <- import("Zea_mays.gtf")
gtf_df <- as.data.frame(gtf)

# Get unique gene coordinates
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

print("4. Building Genotype Matrix from Introgressions...")

# Flatten Teogeno list into GRanges
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

# Find overlaps (which samples have introgression at which gene?)
hits <- findOverlaps(segs_gr, genes_gr)

# Construct matrix (0 = No Introgression, 1 = Introgression)
genotype_mat <- matrix(
  0,
  nrow = nrow(gene_coords),
  ncol = nrow(sample_df),
  dimnames = list(gene_coords$gene_id, sample_df$sample_id)
)

# Fill matrix
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

print("5. Preparing MatrixEQTL Data...")

# Filter non-informative genes
var_mask <- rowSums(genotype_mat) > 0 & rowSums(genotype_mat) < ncol(genotype_mat)
genotype_mat <- genotype_mat[var_mask, , drop = FALSE]

final_genes <- rownames(genotype_mat)
gene_coords_subset <- gene_coords[match(final_genes, gene_coords$gene_id), ]

# Normalize expression (log2 CPM)
lib_size <- colSums(counts)
expr_mat <- as.matrix(counts[rownames(genotype_mat), ])
expr_mat <- log2((t(t(expr_mat) / lib_size)) * 1e6 + 1)

# Create sliced data objects
gene_slice <- SlicedData$new()
gene_slice$CreateFromMatrix(expr_mat)
gene_slice$ResliceCombined(sliceSize = 2000)

snps_slice <- SlicedData$new()
snps_slice$CreateFromMatrix(genotype_mat)
snps_slice$ResliceCombined(sliceSize = 2000)

# Covariates: plate only (remove intercept)
print("   > Processing Covariates (Plate Only)...")
cvrt_mat <- model.matrix(~ plate, data = sample_df)
cvrt_mat <- cvrt_mat[, -1, drop = FALSE]  # Remove intercept
cvrt_mat <- t(cvrt_mat)  # Transpose for MatrixEQTL

print(paste("   > Covariates used:", paste(rownames(cvrt_mat), collapse = ", ")))

cvrt_slice <- SlicedData$new()
cvrt_slice$CreateFromMatrix(cvrt_mat)

# Location files
snpspos <- as.data.frame(gene_coords_subset[, c("gene_id", "chr", "start")])
colnames(snpspos) <- c("snp", "chr", "pos")

genepos <- as.data.frame(gene_coords_subset[, c("gene_id", "chr", "start", "end")])
colnames(genepos) <- c("geneid", "chr", "left", "right")

# ==============================================================================
# 6. RUN MATRIXEQTL
# ==============================================================================

print("6. Running MatrixEQTL...")

# Create output directory if it doesn't exist, check permissions
if (!dir.exists("output")) {
  dir.create("output")
} else {
  # Check if we can write to the directory
  test_file <- "output/.write_test"
  tryCatch({
    writeLines("test", test_file)
    file.remove(test_file)
  }, error = function(e) {
    stop("Output directory exists but is not writable. Check permissions.")
  })
}

me <- Matrix_eQTL_main(
  snps = snps_slice,
  gene = gene_slice,
  cvrt = cvrt_slice,
  output_file_name = NULL,
  pvOutputThreshold = 0,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "output/BZea_cis_eQTL_results_spatial.csv",
  pvOutputThreshold.cis = 1,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 10,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

print("SUCCESS! Results saved to output/BZea_cis_eQTL_results_spatial.csv")

# ==============================================================================
# 7. LOAD AND PREPARE RESULTS
# ==============================================================================

print("7. Loading Results...")

# Load results (tab-delimited)
results <- read.delim(
  "output/BZea_cis_eQTL_results_spatial.csv",
  sep = "\t",
  stringsAsFactors = FALSE
)

# Save Excel-friendly version
write.csv(results, "output/BZea_cis_eQTL_results_spatial_excel.csv", row.names = FALSE)

# Basic statistics
n_sig <- sum(results$FDR < 0.05)
print(paste("Significant hits (FDR < 0.05):", n_sig))
print(paste("Total tests:", nrow(results)))

# Add annotations
results <- results |>
  mutate(
    is_sig = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    log_p = -log10(p.value)
  )

# Save significant results
sig_results <- results |> filter(FDR < 0.05)
write.csv(sig_results, "output/BZea_Significant_eQTLs_spatial.csv", row.names = FALSE)

# Display top 10 hits
top_hits <- results |>
  arrange(p.value) |>
  head(10)

print("Top 10 Significant Genes:")
print(top_hits)

# ==============================================================================
# 8. GENERATE PLOTS
# ==============================================================================

print("8. Generating Plots...")

# ------------------------------------------------------------------------------
# 8a. P-value histogram
# ------------------------------------------------------------------------------

p1 <- ggplot(results, aes(x = p.value)) +
  geom_histogram(bins = 50, fill = "dodgerblue", color = "white") +
  theme_minimal() +
  labs(
    title = "P-value Distribution",
    subtitle = paste("Total Tests:", nrow(results)),
    x = "P-value",
    y = "Frequency"
  )

ggsave("output/plot_pvalue_hist_spatial.png", p1, width = 6, height = 4)

# ------------------------------------------------------------------------------
# 8b. Volcano plot
# ------------------------------------------------------------------------------

p2 <- ggplot(results, aes(x = beta, y = log_p, color = is_sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Introgression Effects",
    subtitle = paste(n_sig, "significant genes (FDR < 0.05)"),
    x = "Effect Size (Beta)",
    y = "-log10(P-value)",
    color = "Status"
  ) +
  theme(legend.position = "top")

ggsave("output/plot_volcano_spatial.png", p2, width = 6, height = 5)

# ------------------------------------------------------------------------------
# 8c. Manhattan plot
# ------------------------------------------------------------------------------

# Prepare data
plot_data <- results |>
  inner_join(snpspos, by = c("SNP" = "snp")) |>
  mutate(
    chr_clean = gsub("chr", "", chr, ignore.case = TRUE), 
    chr_num = as.numeric(chr_clean),
    log_p = -log10(p.value)
  ) |>
  filter(!is.na(chr_num))

# Calculate cumulative positions
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

# Calculate axis labels
axis_set <- plot_data |>
  group_by(chr_num) |>
  summarize(
    center = (as.numeric(max(bp_cum)) + as.numeric(min(bp_cum))) / 2, 
    .groups = "drop"
  ) |>
  arrange(chr_num)

sig_cutoff <- -log10(max(results$p.value[results$FDR < 0.05]))
max_y <- max(plot_data$log_p, na.rm = TRUE)

# Generate plot
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
    title = "Cis-eQTL Manhattan Plot",
    subtitle = paste("Genome-wide view of", nrow(plot_data), "genes"),
    x = "Chromosome",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 10, color = "black")
  )

ggsave("output/ciseqtl_manhattan_spatial.png", p_man, width = 10, height = 5)

# ==============================================================================
# 9. INDIVIDUAL GENE PLOTS
# ==============================================================================

print("9. Generating Individual Gene Plots...")

# ------------------------------------------------------------------------------
# 9a. Top hit with taxa labels
# ------------------------------------------------------------------------------

target_gene <- "Zm00001eb012750"

# Extract data and join metadata
gene_geno <- genotype_mat[target_gene, ]
gene_expr <- expr_mat[target_gene, ]

plot_df <- tibble(
  sample_id = names(gene_geno),
  Genotype = factor(
    ifelse(gene_geno == 1, "Teosinte Introgression", "B73 Background"),
    levels = c("B73 Background", "Teosinte Introgression")
  ),
  Expression = gene_expr
) |>
  left_join(metadata |> select(sample_id, taxa), by = "sample_id") |>
  mutate(
    point_color = ifelse(Genotype == "Teosinte Introgression", taxa, "B73")
  )

# Summary statistics
group_stats <- plot_df |>
  group_by(Genotype) |>
  summarise(
    mean_expr = mean(Expression),
    sd_expr   = sd(Expression),
    .groups = "drop"
  )

# Generate plot
p <- ggplot(plot_df, aes(x = Genotype, y = Expression)) +
  geom_jitter(
    aes(fill = point_color),
    shape = 21,
    size = 3,
    alpha = 0.7,
    color = "grey30",
    position = position_jitter(width = 0.2)
  ) +
  geom_errorbar(
    data = group_stats,
    aes(x = Genotype, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
    width = 0,
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_stats,
    aes(x = Genotype, y = mean_expr),
    shape = 16,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Genotype",
    y = "Normalized Expression (Log2)",
    title = paste("Effect of Introgression on", target_gene),
    subtitle = "Teo points colored by taxa; black = mean ± SD"
  ) +
  scale_x_discrete(
    labels = c("B73 Background" = "B73", "Teosinte Introgression" = "Teo")
  ) +
  scale_fill_manual(
    values = c(
      "B73"  = "#03bec4",
      "Bals" = "#f364e2",
      "Zdip" = "#f8756d",
      "Hueh" = "#b69d00",
      "Zlux" = "#00b837",
      "Nobo" = "#609bfe",
      "Mesa" = "#609bfe",
      "Chal" = "#609bfe",
      "Dura" = "#609bfe",
      "Nabo" = "#609bfe"
    ),
    name = "Teosinte taxa"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

ggsave("output/eQTL_Zm00001eb012750_taxa_spatial.png", p, width = 5, height = 5)

# ------------------------------------------------------------------------------
# 9b. Top 10 hits panel
# ------------------------------------------------------------------------------

top10_genes <- c(
  "Zm00001eb012750",
  "Zm00001eb375600",
  "Zm00001eb401220",
  "Zm00001eb157760",
  "Zm00001eb039390",
  "Zm00001eb015830",
  "Zm00001eb077910",
  "Zm00001eb108910",
  "Zm00001eb013080",
  "Zm00001eb399300"
)

# Build long dataframe
plot_df_all <- map_dfr(top10_genes, function(gene) {
  gene_geno <- genotype_mat[gene, ]
  gene_expr <- expr_mat[gene, ]
  
  tibble(
    Gene = gene,
    Sample = names(gene_geno),
    Genotype = factor(
      ifelse(gene_geno == 1, "Teosinte Introgression", "B73 Background"),
      levels = c("B73 Background", "Teosinte Introgression")
    ),
    Expression = gene_expr
  )
})

# Summary stats
group_stats_all <- plot_df_all |>
  group_by(Gene, Genotype) |>
  summarise(
    mean_expr = mean(Expression),
    sd_expr   = sd(Expression),
    .groups = "drop"
  )

# Paneled plot
p <- ggplot(plot_df_all, aes(x = Genotype, y = Expression, fill = Genotype)) +
  geom_jitter(
    shape = 21,
    size = 2,
    alpha = 0.6,
    color = "white",
    position = position_jitter(width = 0.2)
  ) +
  geom_errorbar(
    data = group_stats_all,
    aes(x = Genotype, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
    width = 0,
    size = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_stats_all,
    aes(x = Genotype, y = mean_expr),
    shape = 16,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  labs(
    x = "Genotype",
    y = "Normalized Expression (Log2)",
    title = "Top 10 cis-eQTL Hits",
    subtitle = "Black = Mean ± SD"
  ) +
  scale_fill_manual(
    values = c(
      "B73 Background" = "grey70",
      "Teosinte Introgression" = "forestgreen"
    ),
    labels = c("B73 Background" = "B73", "Teosinte Introgression" = "Teo")
  ) +
  scale_x_discrete(
    labels = c("B73 Background" = "B73", "Teosinte Introgression" = "Teo")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave("output/eQTL_top10_hits_panel_spatial.png", p, width = 12, height = 8)

# ==============================================================================
# 10. TRANS-EQTL FOR A SINGLE GENE (EXCLUDING CIS)
# ==============================================================================
# trans-eqtl scan for: Zm00001eb372490
# coordinates: chr9 6311753 6313844
# note: we exclude markers in a cis window around the gene, then test all remaining markers genome-wide

print("10. Running trans-eQTL (excluding cis) for Zm00001eb372490")

target_gene_trans  <- "Zm00001eb372490"
target_chr_trans   <- "chr9"
target_start_trans <- 6311753
target_end_trans   <- 6313844

# ---- cis exclusion window (bp) ----
# adjust if you want a wider/narrower exclusion around the gene
cis_exclusion_bp <- 1e6

cis_left  <- target_start_trans - cis_exclusion_bp
cis_right <- target_end_trans   + cis_exclusion_bp

# ---- build marker position table for all genotype markers currently in genotype_mat ----
marker_coords <- gene_coords |>
  filter(gene_id %in% rownames(genotype_mat)) |>
  transmute(
    snp = gene_id,
    chr = as.character(chr),
    pos = as.integer(start)
  )

# ---- drop cis markers near the target gene on chr9 ----
cis_marker_mask <- marker_coords$chr == target_chr_trans &
  marker_coords$pos >= cis_left &
  marker_coords$pos <= cis_right

snps_keep <- marker_coords$snp[!cis_marker_mask]

# also drop the target gene marker itself if it exists as a marker row
snps_keep <- setdiff(snps_keep, target_gene_trans)

genotype_mat_trans <- genotype_mat[snps_keep, , drop = FALSE]
snpspos_trans <- marker_coords |>
  filter(snp %in% snps_keep)

print(paste0(
  "   > kept markers for trans test: ",
  nrow(genotype_mat_trans),
  " (excluded cis window +/- ",
  cis_exclusion_bp,
  " bp)"
))

# ---- expression slice: target gene only ----
expr_target <- expr_mat[target_gene_trans, , drop = FALSE]

gene_slice_trans <- SlicedData$new()
gene_slice_trans$CreateFromMatrix(expr_target)
gene_slice_trans$ResliceCombined(sliceSize = 2000)

snps_slice_trans <- SlicedData$new()
snps_slice_trans$CreateFromMatrix(genotype_mat_trans)
snps_slice_trans$ResliceCombined(sliceSize = 2000)

# ---- gene position table (required if you want MatrixEQTL to reason about cis; here it's mostly bookkeeping) ----
genepos_trans <- tibble(
  geneid = target_gene_trans,
  chr    = target_chr_trans,
  left   = target_start_trans,
  right  = target_end_trans
)

# ---- run MatrixEQTL for trans only (cis already excluded by filtering markers) ----
me_trans <- Matrix_eQTL_main(
  snps = snps_slice_trans,
  gene = gene_slice_trans,
  cvrt = cvrt_slice,                       # re-use your plate covariates
  output_file_name = "output/BZea_trans_eQTL_Zm00001eb372490_excl_cis.tsv",
  pvOutputThreshold = 1,                   # keep everything; filter later if you want
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = NULL,             # do not write a cis file for this run
  pvOutputThreshold.cis = 0,
  snpspos = snpspos_trans,
  genepos = genepos_trans,
  cisDist = cis_exclusion_bp,              # not relied on for exclusion (we already filtered), but consistent
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

print("   > trans-eQTL run complete!")

# ---- load + save excel-friendly copy ----
trans_results <- read.delim(
  "output/BZea_trans_eQTL_Zm00001eb372490_excl_cis.tsv",
  sep = "\t",
  stringsAsFactors = FALSE
)

write.csv(
  trans_results,
  "output/BZea_trans_eQTL_Zm00001eb372490_excl_cis_excel.csv",
  row.names = FALSE
)

# ---- quick summary ----
trans_results <- trans_results |>
  mutate(
    is_sig = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    log_p  = -log10(p.value)
  )

n_sig_trans <- sum(trans_results$FDR < 0.05)
print(paste("   > significant trans hits (FDR < 0.05):", n_sig_trans))
print(paste("   > total trans tests:", nrow(trans_results)))

top_hits_trans <- trans_results |>
  arrange(p.value) |>
  head(10)

print("Top 10 trans hits:")
print(top_hits_trans)

# ==============================================================================
# 10b. TRANS MANHATTAN PLOT (Zm00001eb372490, cis excluded)
# ==============================================================================

print("10b. Plotting trans Manhattan...")

trans_plot <- trans_results |>
  inner_join(snpspos_trans, by = c("SNP" = "snp")) |>
  mutate(
    chr_clean = gsub("chr", "", chr, ignore.case = TRUE),
    chr_num   = as.numeric(chr_clean),
    pos       = as.numeric(pos),
    log_p     = -log10(p.value)
  ) |>
  filter(!is.na(chr_num), !is.na(pos))

# cumulative x positions
trans_cum <- trans_plot |>
  group_by(chr_num) |>
  summarise(max_bp = max(pos), .groups = "drop") |>
  arrange(chr_num) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr_num, bp_add)

trans_plot <- trans_plot |>
  inner_join(trans_cum, by = "chr_num") |>
  mutate(
    bp_cum = pos + bp_add,
    color_group = case_when(
      FDR < 0.05 & (chr_num %% 2 == 1) ~ "Sig_Odd",
      FDR < 0.05 & (chr_num %% 2 == 0) ~ "Sig_Even",
      (chr_num %% 2 == 1) ~ "Base_Odd",
      TRUE ~ "Base_Even"
    )
  )

axis_set_trans <- trans_plot |>
  group_by(chr_num) |>
  summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop") |>
  arrange(chr_num)

sig_cutoff_trans <- if (any(trans_results$FDR < 0.05)) {
  -log10(max(trans_results$p.value[trans_results$FDR < 0.05]))
} else {
  NA_real_
}

max_y_trans <- max(trans_plot$log_p, na.rm = TRUE)

p_man_trans <- ggplot(trans_plot, aes(x = bp_cum, y = log_p, color = color_group)) +
  geom_point(alpha = 0.75, size = 1.3) +
  scale_color_manual(values = c(
    "Base_Odd"  = "grey45",
    "Base_Even" = "grey70",
    "Sig_Odd"   = "firebrick",
    "Sig_Even"  = "darkred"
  )) +
  { if (!is.na(sig_cutoff_trans)) geom_hline(yintercept = sig_cutoff_trans, color = "red", linetype = "dashed") } +
  scale_x_continuous(
    breaks = axis_set_trans$center,
    labels = as.character(axis_set_trans$chr_num)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max_y_trans * 1.1)
  ) +
  labs(
    title = "Trans-eQTL Manhattan Plot",
    subtitle = paste0(
      "Target gene: Zm00001eb372490; cis excluded +/- ",
      format(cis_exclusion_bp, scientific = FALSE),
      " bp; tests = ", nrow(trans_plot)
    ),
    x = "Chromosome",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 10, color = "black")
  )

p_man_trans

ggsave("output/trans_manhattan_tcptf9_excl_cis.png", p_man_trans, width = 10, height = 5)

# ==============================================================================
# 10c. TRANS-EQTL + MANHATTAN FOR MULTIPLE TARGET GENES (EXCLUDING CIS)
# ==============================================================================

print("10c. Running trans-eQTL + Manhattan for multiple targets...")

# cis exclusion window (bp)
cis_exclusion_bp <- 1e6

# marker position table for all genotype markers
marker_coords <- gene_coords |>
  filter(gene_id %in% rownames(genotype_mat)) |>
  transmute(
    snp = gene_id,
    chr = as.character(chr),
    pos = as.integer(start)
  )

run_trans_and_manhattan <- function(gene_name, gene_id, gene_chr, gene_start, gene_end) {
  
  print(paste0("   > ", gene_name, " (", gene_id, "): trans scan starting..."))
  
  cis_left  <- gene_start - cis_exclusion_bp
  cis_right <- gene_end   + cis_exclusion_bp
  
  # exclude cis markers around this gene
  cis_marker_mask <- marker_coords$chr == gene_chr &
    marker_coords$pos >= cis_left &
    marker_coords$pos <= cis_right
  
  snps_keep <- marker_coords$snp[!cis_marker_mask]
  snps_keep <- setdiff(snps_keep, gene_id)
  
  genotype_mat_trans <- genotype_mat[snps_keep, , drop = FALSE]
  snpspos_trans <- marker_coords |>
    filter(snp %in% snps_keep)
  
  print(paste0(
    "      kept markers: ", nrow(genotype_mat_trans),
    " (excluded cis +/- ", format(cis_exclusion_bp, scientific = FALSE), " bp)"
  ))
  
  # expression: target gene only
  expr_target <- expr_mat[gene_id, , drop = FALSE]
  
  gene_slice_trans <- SlicedData$new()
  gene_slice_trans$CreateFromMatrix(expr_target)
  gene_slice_trans$ResliceCombined(sliceSize = 2000)
  
  snps_slice_trans <- SlicedData$new()
  snps_slice_trans$CreateFromMatrix(genotype_mat_trans)
  snps_slice_trans$ResliceCombined(sliceSize = 2000)
  
  genepos_trans <- tibble(
    geneid = gene_id,
    chr    = gene_chr,
    left   = gene_start,
    right  = gene_end
  )
  
  out_tsv <- paste0("output/BZea_trans_eQTL_", gene_name, "_excl_cis.tsv")
  out_csv <- paste0("output/BZea_trans_eQTL_", gene_name, "_excl_cis_excel.csv")
  out_png <- paste0("output/trans_manhattan_", gene_name, "_excl_cis.png")
  
  me_trans <- Matrix_eQTL_main(
    snps = snps_slice_trans,
    gene = gene_slice_trans,
    cvrt = cvrt_slice,
    output_file_name = out_tsv,
    pvOutputThreshold = 1,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 0,
    snpspos = snpspos_trans,
    genepos = genepos_trans,
    cisDist = cis_exclusion_bp,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  trans_results <- read.delim(out_tsv, sep = "\t", stringsAsFactors = FALSE)
  
  write.csv(trans_results, out_csv, row.names = FALSE)
  
  trans_results <- trans_results |>
    mutate(
      is_sig = ifelse(FDR < 0.05, "Significant", "Not Significant"),
      log_p  = -log10(p.value)
    )
  
  n_sig_trans <- sum(trans_results$FDR < 0.05)
  print(paste0("      significant (FDR < 0.05): ", n_sig_trans, " / ", nrow(trans_results)))
  
  # ---- Manhattan prep ----
  trans_plot <- trans_results |>
    inner_join(snpspos_trans, by = c("SNP" = "snp")) |>
    mutate(
      chr_clean = gsub("chr", "", chr, ignore.case = TRUE),
      chr_num   = as.numeric(chr_clean),
      pos       = as.numeric(pos),
      log_p     = -log10(p.value)
    ) |>
    filter(!is.na(chr_num), !is.na(pos))
  
  trans_cum <- trans_plot |>
    group_by(chr_num) |>
    summarise(max_bp = max(pos), .groups = "drop") |>
    arrange(chr_num) |>
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
    select(chr_num, bp_add)
  
  trans_plot <- trans_plot |>
    inner_join(trans_cum, by = "chr_num") |>
    mutate(
      bp_cum = pos + bp_add,
      color_group = case_when(
        FDR < 0.05 & (chr_num %% 2 == 1) ~ "Sig_Odd",
        FDR < 0.05 & (chr_num %% 2 == 0) ~ "Sig_Even",
        (chr_num %% 2 == 1) ~ "Base_Odd",
        TRUE ~ "Base_Even"
      )
    )
  
  axis_set_trans <- trans_plot |>
    group_by(chr_num) |>
    summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop") |>
    arrange(chr_num)
  
  sig_cutoff_trans <- if (any(trans_results$FDR < 0.05)) {
    -log10(max(trans_results$p.value[trans_results$FDR < 0.05]))
  } else {
    NA_real_
  }
  
  max_y_trans <- max(trans_plot$log_p, na.rm = TRUE)
  
  p_man_trans <- ggplot(trans_plot, aes(x = bp_cum, y = log_p, color = color_group)) +
    geom_point(alpha = 0.75, size = 1.3) +
    scale_color_manual(values = c(
      "Base_Odd"  = "grey45",
      "Base_Even" = "grey70",
      "Sig_Odd"   = "firebrick",
      "Sig_Even"  = "darkred"
    )) +
    { if (!is.na(sig_cutoff_trans)) geom_hline(yintercept = sig_cutoff_trans, color = "red", linetype = "dashed") } +
    scale_x_continuous(
      breaks = axis_set_trans$center,
      labels = as.character(axis_set_trans$chr_num)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max_y_trans * 1.1)
    ) +
    labs(
      title = paste0("Trans-eQTL Manhattan: ", gene_name),
      subtitle = paste0(
        gene_id, " (", gene_chr, ":", gene_start, "-", gene_end, "); cis excluded +/- ",
        format(cis_exclusion_bp, scientific = FALSE), " bp; tests = ", nrow(trans_plot)
      ),
      x = "Chromosome",
      y = "-log10(P-value)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 10, color = "black")
    )
  
  p_man_trans
  ggsave(out_png, p_man_trans, width = 10, height = 5)
  
  print(paste0("      saved: ", out_png))
  
  return(invisible(list(results = trans_results, plot = p_man_trans)))
}

# targets
targets <- tribble(
  ~gene_name, ~gene_id,           ~chr,   ~start,    ~end,
  "hpc1",     "Zm00001eb121780",  "chr3",  8496140,   8501237,
  "nlp1",     "Zm00001eb231720",  "chr5",  78778520,  78781733,
  "nrg11",    "Zm00001eb206940",  "chr4",  244043231, 244045806
)

# run all
trans_runs <- pmap(
  targets,
  ~run_trans_and_manhattan(..1, ..2, ..3, ..4, ..5)
)

print("10c. Done: trans Manhattan plots saved for hpc1, nlp1, nrg11.")

