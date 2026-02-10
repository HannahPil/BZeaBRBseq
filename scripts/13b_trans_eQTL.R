# ==============================================================================
# BZEA RNA: TRANS-EQTL PIPELINE (CHROMOSOME-BY-CHROMOSOME)
# ==============================================================================
# purpose:
#   run trans-eQTLs using gene-level introgression genotypes (0/1) as predictors
#   and log2 CPM expression as outcomes, controlling for plate.
#
# strategy (memory-safe):
#   for each source chromosome:
#     - snps = introgression status for genes on that chromosome (rows)
#     - genes = expression for all genes with genotype variation (rows)
#     - run MatrixEQTL with pvOutputThreshold to write only small p-values
#     - remove cis-by-distance results from the saved file
#
# outputs:
#   output/trans_by_source_chr/
#     BZea_trans_sources_chr1_pv1e-04.tsv (+ .csv)
#     ...
#
# notes:
# - trans definition here: NOT diagonal and NOT within cisDist_bp of the target gene interval
# - requires: Zea_mays.gtf, results_list_new_name.rds, Zea_mays_counts.txt, metadata.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MatrixEQTL)
  library(rtracklayer)
  library(GenomicRanges)
})

# ==============================================================================
# 0. SETTINGS
# ==============================================================================
out_dir <- file.path("output", "trans_by_source_chr")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pv_trans   <- 1e-4   # p-value threshold for saving trans hits (controls file size)
cisDist_bp <- 1e6    # cis exclusion window around target gene interval (bp)
slice_size <- 2000

cat("settings:\n")
cat("  out_dir     = ", out_dir, "\n", sep = "")
cat("  pv_trans    = ", pv_trans, "\n", sep = "")
cat("  cisDist_bp  = ", format(cisDist_bp, scientific = FALSE), "\n", sep = "")
cat("  slice_size  = ", slice_size, "\n", sep = "")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("\n1. Loading data...\n")

teogeno  <- readRDS("results_list_new_name.rds")
counts   <- read.delim("Zea_mays_counts.txt", check.names = FALSE, row.names = 1)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

teogeno <- teogeno[!duplicated(names(teogeno))]

# ==============================================================================
# 2. ALIGN SAMPLES
# ==============================================================================
cat("\n2. Aligning samples...\n")

sample_df <- metadata |>
  filter(sample_id %in% colnames(counts)) |>
  transmute(
    sample_id,
    plate = factor(plate),
    genotype_teogeno_key = if_else(genotype == "B73", "B73.B", paste0(genotype, ".B"))
  ) |>
  filter(genotype_teogeno_key %in% names(teogeno)) |>
  arrange(sample_id)

counts <- counts[, as.character(sample_df$sample_id), drop = FALSE]

cat("  samples kept: ", nrow(sample_df), "\n", sep = "")

# ==============================================================================
# 3. GENE COORDINATES (GTF)
# ==============================================================================
cat("\n3. Importing GTF + building gene coordinates...\n")

gtf <- import("Zea_mays.gtf")
gtf_df <- as.data.frame(gtf)

gene_coords <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(
    chr   = as.character(seqnames[1]),
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

cat("  genes with coords + counts: ", nrow(gene_coords), "\n", sep = "")

# ==============================================================================
# 4. BUILD GENE-LEVEL INTROGRESSION GENOTYPE MATRIX (0/1 PER GENE x SAMPLE)
# ==============================================================================
cat("\n4. Building gene-level introgression matrix...\n")

seg_df <- imap_dfr(teogeno, ~{
  .x |>
    filter(V4 == "Introgression") |>
    transmute(
      key   = .y,
      chr   = V1,
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

geneGT <- matrix(
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
    geneGT[g_id, as.character(s_ids)] <- 1
  }
}

cat("  introgressed gene markers (any 1s): ", sum(rowSums(geneGT) > 0), "\n", sep = "")

# ==============================================================================
# 5. NORMALIZE EXPRESSION (LOG2 CPM)
# ==============================================================================
cat("\n5. Normalizing expression (log2 CPM)...\n")

lib_size <- colSums(counts)
expr_mat <- log2(t(t(as.matrix(counts)) / lib_size * 1e6) + 1)

# ==============================================================================
# 6. SUBSET TO GENES WITH GENOTYPE VARIATION
# ==============================================================================
cat("\n6. Subsetting to genes with genotype variation...\n")

var_genes <- rownames(geneGT)[rowSums(geneGT) > 0 & rowSums(geneGT) < ncol(geneGT)]
var_genes <- intersect(var_genes, rownames(expr_mat))

cat("  genes with genotype variation: ", length(var_genes), "\n", sep = "")

geneGT_var <- geneGT[var_genes, , drop = FALSE]
expr_var   <- expr_mat[var_genes, , drop = FALSE]
gene_coords_var <- gene_coords |>
  filter(gene_id %in% var_genes)

# ==============================================================================
# 7. COVARIATES (PLATE)
# ==============================================================================
cat("\n7. Building covariates (plate)...\n")

cvrt_mat <- model.matrix(~ plate, data = sample_df)
cvrt_mat <- cvrt_mat[, -1, drop = FALSE]
cvrt_mat <- t(cvrt_mat)

cvrt_slice <- SlicedData$new()
cvrt_slice$CreateFromMatrix(cvrt_mat)

cat("  covariates: ", paste(rownames(cvrt_mat), collapse = ", "), "\n", sep = "")

# ==============================================================================
# 8. PRECOMPUTE TARGET GENE POSITION TABLE (FOR CIS EXCLUSION)
# ==============================================================================
cat("\n8. Building target gene position table...\n")

genepos_all <- tibble(geneid = rownames(expr_var)) |>
  left_join(
    gene_coords_var |> select(gene_id, chr, start, end),
    by = c("geneid" = "gene_id")
  ) |>
  transmute(
    geneid = geneid,
    chr    = as.character(chr),
    left   = as.integer(start),
    right  = as.integer(end)
  ) |>
  as.data.frame()

stopifnot(!any(is.na(genepos_all$chr)), !any(is.na(genepos_all$left)), !any(is.na(genepos_all$right)))
genepos_all <- genepos_all[match(rownames(expr_var), genepos_all$geneid), ]

# ==============================================================================
# 9. RUN TRANS EQTLS CHR-BY-CHR (SOURCE CHR)
# ==============================================================================
cat("\n9. Running trans-eQTLs by source chromosome...\n")

chr_levels <- gene_coords_var |>
  distinct(chr) |>
  mutate(chr_num = as.numeric(gsub("chr", "", chr, ignore.case = TRUE))) |>
  arrange(chr_num) |>
  pull(chr)

run_trans_chr <- function(chr_run) {
  
  cat("\n  -> sources on ", chr_run, "\n", sep = "")
  
  snps_chr <- gene_coords_var |>
    filter(chr == chr_run) |>
    pull(gene_id)
  
  if(length(snps_chr) == 0) {
    cat("     skipping: no source genes\n")
    return(invisible(NULL))
  }
  
  geneGT_chr <- geneGT_var[snps_chr, , drop = FALSE]
  
  # snps slice
  snps_slice <- SlicedData$new()
  snps_slice$CreateFromMatrix(as.matrix(geneGT_chr))
  snps_slice$ResliceCombined(sliceSize = slice_size)
  
  # genes slice (targets = all var genes)
  gene_slice <- SlicedData$new()
  gene_slice$CreateFromMatrix(as.matrix(expr_var))
  gene_slice$ResliceCombined(sliceSize = slice_size)
  
  # build snpspos (must match rownames of geneGT_chr)
  snpspos_chr <- tibble(SNP = rownames(geneGT_chr)) |>
    left_join(
      gene_coords_var |> select(gene_id, chr, start),
      by = c("SNP" = "gene_id")
    ) |>
    transmute(
      SNP = SNP,
      SNP_chr = as.character(chr),
      SNP_pos = as.integer(start)
    )
  
  stopifnot(!any(is.na(snpspos_chr$SNP_chr)), !any(is.na(snpspos_chr$SNP_pos)))
  snpspos_chr <- snpspos_chr[match(rownames(geneGT_chr), snpspos_chr$SNP), ]
  
  # write snpspos for MatrixEQTL (it requires cols named snp/chr/pos)
  snpspos_for_meqtl <- snpspos_chr |>
    transmute(snp = SNP, chr = SNP_chr, pos = SNP_pos) |>
    as.data.frame()
  
  out_trans_raw <- file.path(out_dir, paste0("BZea_trans_sources_", chr_run, "_pv", pv_trans, "_RAW.tsv"))
  out_trans     <- file.path(out_dir, paste0("BZea_trans_sources_", chr_run, "_pv", pv_trans, ".tsv"))
  
  Matrix_eQTL_main(
    snps = snps_slice,
    gene = gene_slice,
    cvrt = cvrt_slice,
    output_file_name = out_trans_raw,
    pvOutputThreshold = pv_trans,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 0,
    snpspos = snpspos_for_meqtl,
    genepos = genepos_all,
    cisDist = cisDist_bp,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  # load, drop diagonal
  trans_res <- read.delim(out_trans_raw, sep = "\t", stringsAsFactors = FALSE) |>
    filter(SNP != gene)
  
  # build genepos join table with final names up front (no renaming)
  genepos_join <- as_tibble(genepos_all) |>
    transmute(
      gene = geneid,
      gene_chr = as.character(chr),
      gene_left = as.integer(left),
      gene_right = as.integer(right)
    )
  
  # drop cis-by-distance
  trans_res <- trans_res |>
    left_join(snpspos_chr, by = "SNP") |>
    left_join(genepos_join, by = "gene") |>
    filter(!(SNP_chr == gene_chr &
               SNP_pos >= (gene_left - cisDist_bp) &
               SNP_pos <= (gene_right + cisDist_bp))) |>
    select(-SNP_chr, -SNP_pos, -gene_chr, -gene_left, -gene_right)
  
  # write cleaned outputs
  write.table(trans_res, out_trans, sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(trans_res, sub("\\.tsv$", ".csv", out_trans), row.names = FALSE)
  
  file.remove(out_trans_raw)
  
  cat("     saved: ", out_trans, "\n", sep = "")
  cat("     n kept: ", nrow(trans_res), "\n", sep = "")
  cat("     sig (FDR < 0.05): ", sum(trans_res$FDR < 0.05), "\n", sep = "")
  
  invisible(NULL)
}


for(chr_run in chr_levels) {
  run_trans_chr(chr_run)
}

cat("\nDONE. trans results in: ", out_dir, "\n", sep = "")

# ==============================================================================
# PLOT TRANS MANHATTAN FOR ALL SOURCE-CHR RUNS
# ==============================================================================
# reads per-chromosome trans result files in output/trans_by_source_chr/
# and makes one Manhattan per file, saving PNGs into the same folder.
#
# x-axis = genomic position of SOURCE gene (the "SNP" gene)
# y-axis = -log10(p.value)
#
# assumptions:
# - you ran the trans loop and produced files like:
#     output/trans_by_source_chr/BZea_trans_sources_chr1_pv1e-04.tsv
# - you still have gene_coords_var in memory (or we rebuild it below)
# ==============================================================================

cat("Plotting trans Manhattan plots for all chromosomes...\n")

plot_dir <- file.path("output", "trans_by_source_chr")
stopifnot(dir.exists(plot_dir))

# if gene_coords_var is not in memory, rebuild this minimal table from gene_coords
# (comment this out if you already have gene_coords_var)
# gene_coords_var <- gene_coords |> filter(gene_id %in% rownames(expr_var))

tsv_files <- list.files(
  plot_dir,
  pattern = "^BZea_trans_sources_chr[0-9]+_pv.*\\.tsv$",
  full.names = TRUE
)

cat("  found ", length(tsv_files), " trans tsv files\n", sep = "")

# helper to make cumulative x positions
make_cumpos <- function(df) {
  df_cum <- df |>
    group_by(chr_num) |>
    summarise(max_bp = max(pos, na.rm = TRUE), .groups = "drop") |>
    arrange(chr_num) |>
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
    select(chr_num, bp_add)
  
  df |>
    inner_join(df_cum, by = "chr_num") |>
    mutate(bp_cum = pos + bp_add)
}

for(f in tsv_files) {
  
  chr_run <- sub(".*(chr[0-9]+)_pv.*", "\\1", basename(f))
  cat("  -> ", chr_run, "\n", sep = "")
  
  trans_res <- read.delim(f, sep = "\t", stringsAsFactors = FALSE)
  
  if(nrow(trans_res) == 0) {
    cat("     (skipping plot: no rows)\n")
    next
  }
  
  # join SOURCE (SNP) positions from gene_coords_var
  plot_data <- trans_res |>
    left_join(
      gene_coords_var |> select(gene_id, chr, start),
      by = c("SNP" = "gene_id")
    ) |>
    mutate(
      chr_clean = gsub("chr", "", chr, ignore.case = TRUE),
      chr_num   = suppressWarnings(as.numeric(chr_clean)),
      pos       = as.numeric(start),
      log_p     = -log10(p.value)
    ) |>
    filter(!is.na(chr_num), !is.na(pos), !is.na(log_p))
  
  if(nrow(plot_data) == 0) {
    cat("     (skipping plot: no rows after coord join)\n")
    next
  }
  
  plot_data <- make_cumpos(plot_data)
  
  axis_set <- plot_data |>
    group_by(chr_num) |>
    summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop") |>
    arrange(chr_num)
  
  max_y <- max(plot_data$log_p, na.rm = TRUE)
  
  sig_cutoff <- if(any(plot_data$FDR < 0.05, na.rm = TRUE)) {
    -log10(max(plot_data$p.value[plot_data$FDR < 0.05], na.rm = TRUE))
  } else {
    NA_real_
  }
  
  plot_data <- plot_data |>
    mutate(
      color_group = case_when(
        FDR < 0.05 & (chr_num %% 2 == 1) ~ "Sig_Odd",
        FDR < 0.05 & (chr_num %% 2 == 0) ~ "Sig_Even",
        (chr_num %% 2 == 1) ~ "Base_Odd",
        TRUE ~ "Base_Even"
      )
    )
  
  p <- ggplot(plot_data, aes(x = bp_cum, y = log_p, color = color_group)) +
    geom_point(alpha = 0.75, size = 1.2) +
    scale_color_manual(values = c(
      "Base_Odd"  = "grey45",
      "Base_Even" = "grey70",
      "Sig_Odd"   = "firebrick",
      "Sig_Even"  = "darkred"
    )) +
    { if (!is.na(sig_cutoff)) geom_hline(yintercept = sig_cutoff, color = "red", linetype = "dashed") } +
    scale_x_continuous(
      breaks = axis_set$center,
      labels = as.character(axis_set$chr_num)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y * 1.1)) +
    labs(
      title = paste0("Trans-eQTL Manhattan (sources = ", chr_run, ")"),
      subtitle = paste0("shown: p <= ", pv_trans, "; n = ", nrow(plot_data)),
      x = "Chromosome (source gene position)",
      y = "-log10(P-value)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 10, color = "black")
    )
  
  out_png <- file.path(plot_dir, paste0("plot_trans_manhattan_sources_", chr_run, "_pv", pv_trans, ".png"))
  ggsave(out_png, p, width = 10, height = 5)
  
  # print to viewer too
  print(p)
}

cat("DONE. plots saved in: ", plot_dir, "\n", sep = "")

# ==============================================================================
# 10. COMBINE ALL TRANS RESULTS + ONE GENOME-WIDE TRANS MANHATTAN
# ==============================================================================
cat("10. Combining trans results + plotting combined Manhattan...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

plot_dir <- file.path("output", "trans_by_source_chr")

tsv_files <- list.files(
  plot_dir,
  pattern = "^BZea_trans_sources_chr[0-9]+_pv.*\\.tsv$",
  full.names = TRUE
)

# read + bind all per-chr trans results
trans_all <- map_dfr(tsv_files, function(f) {
  chr_run <- sub(".*(chr[0-9]+)_pv.*", "\\1", basename(f))
  read.delim(f, sep = "\t", stringsAsFactors = FALSE) |>
    mutate(source_chr_run = chr_run)
})

cat("  total rows combined: ", nrow(trans_all), "\n", sep = "")

# write combined tables
out_tsv <- file.path(plot_dir, "BZea_trans_ALL_sources_combined.tsv")
out_csv <- file.path(plot_dir, "BZea_trans_ALL_sources_combined.csv")

write.table(trans_all, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(trans_all, out_csv, row.names = FALSE)

cat("  saved combined tables:\n")
cat("   - ", out_tsv, "\n", sep = "")
cat("   - ", out_csv, "\n", sep = "")

# ------------------------------------------------------------------------------
# combined trans Manhattan (x = source gene position genome-wide)
# ------------------------------------------------------------------------------

plot_data <- trans_all |>
  left_join(
    gene_coords_var |> select(gene_id, chr, start),
    by = c("SNP" = "gene_id")
  ) |>
  mutate(
    chr_clean = gsub("chr", "", chr, ignore.case = TRUE),
    chr_num   = suppressWarnings(as.numeric(chr_clean)),
    pos       = as.numeric(start),
    log_p     = -log10(p.value)
  ) |>
  filter(!is.na(chr_num), !is.na(pos), !is.na(log_p))

# cumulative positions across genome
data_cum <- plot_data |>
  group_by(chr_num) |>
  summarise(max_bp = max(pos), .groups = "drop") |>
  arrange(chr_num) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr_num, bp_add)

plot_data <- plot_data |>
  inner_join(data_cum, by = "chr_num") |>
  mutate(
    bp_cum = pos + bp_add,
    color_group = case_when(
      FDR < 0.05 & (chr_num %% 2 == 1) ~ "Sig_Odd",
      FDR < 0.05 & (chr_num %% 2 == 0) ~ "Sig_Even",
      (chr_num %% 2 == 1) ~ "Base_Odd",
      TRUE ~ "Base_Even"
    )
  )

axis_set <- plot_data |>
  group_by(chr_num) |>
  summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop") |>
  arrange(chr_num)

max_y <- max(plot_data$log_p, na.rm = TRUE)

sig_cutoff <- if (any(plot_data$FDR < 0.05, na.rm = TRUE)) {
  -log10(max(plot_data$p.value[plot_data$FDR < 0.05], na.rm = TRUE))
} else {
  NA_real_
}

p_trans_all <- ggplot(plot_data, aes(x = bp_cum, y = log_p, color = color_group)) +
  geom_point(alpha = 0.75, size = 1.1) +
  scale_color_manual(values = c(
    "Base_Odd"  = "grey45",
    "Base_Even" = "grey70",
    "Sig_Odd"   = "firebrick",
    "Sig_Even"  = "darkred"
  )) +
  { if (!is.na(sig_cutoff)) geom_hline(yintercept = sig_cutoff, color = "red", linetype = "dashed") } +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = as.character(axis_set$chr_num)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max_y * 1.1)
  ) +
  labs(
    title = "Trans-eQTL Manhattan (all source chromosomes combined)",
    subtitle = paste0(
      "shown are saved hits only (p-value thresholded during scan); n = ",
      nrow(plot_data)
    ),
    x = "Chromosome (source gene position)",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 10, color = "black")
  )

p_trans_all

out_png <- file.path(plot_dir, "plot_trans_manhattan_ALL_sources_combined.png")
ggsave(out_png, p_trans_all, width = 10, height = 5)

cat("  saved plot:\n")
cat("   - ", out_png, "\n", sep = "")

