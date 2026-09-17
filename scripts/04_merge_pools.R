#!/usr/bin/env Rscript

# ==============================================================================
# 04 -- Merge per-pool STARsolo matrices into unified count tables
#
# STARsolo produces one Solo.out/Gene/raw/ directory per pool with:
#   umiDedup-1MM_Directional.mtx    UMI-collapsed counts (canonical)
#   umiDedup-NoDedup.mtx            raw read counts
#   features.tsv                    gene IDs (same across pools; sanity-check)
#   barcodes.tsv                    barcodes actually seen in this pool
#
# This script:
#   - reads all four pools' matrices
#   - maps barcode -> sample_id using the per-pool barcode maps from
#     02_prepare_barcodes.R
#   - column-binds pools into single gene x sample matrices
#   - writes two files in the same shape as the legacy Zea_mays_counts.txt:
#         data/processed/Zea_mays_counts.txt          (UMI-collapsed, canonical)
#         data/processed/Zea_mays_counts_raw.txt      (no dedup, comparison)
#
# Run LOCALLY after copying per-pool Solo.out/Gene/raw/ trees down from HPC:
#   data/starsolo/pool_1/Solo.out/Gene/raw/
#   ... (same for pools 2-4)
# (or, following our git-first convention, have step 03 write directly into
#  hannah/BZeaBRBseq/data/starsolo/ and pull via git.)
# ==============================================================================

suppressPackageStartupMessages(library(Matrix))

data_dir  <- "data"
# STARSOLO_ROOT env var points at the .mtx-holding tree on HPC. Default
# to the local mirror under data/starsolo/ for local runs.
solo_root <- Sys.getenv("STARSOLO_ROOT",
                        unset = file.path(data_dir, "starsolo"))
map_dir   <- file.path(data_dir, "starsolo")     # barcode maps always in repo
out_dir   <- file.path(data_dir, "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("STARsolo output root: ", solo_root, "\n", sep = "")

# --- read per-pool matrix + barcode map -----------------------------------
read_pool <- function(p, matrix_name) {
  raw_dir <- file.path(solo_root, sprintf("pool_%d", p),
                       "Solo.out", "Gene", "raw")
  mtx_path <- file.path(raw_dir, matrix_name)
  bc_path  <- file.path(raw_dir, "barcodes.tsv")
  ft_path  <- file.path(raw_dir, "features.tsv")
  map_path <- file.path(map_dir, sprintf("pool_%d_barcode_map.tsv", p))
  for (f in c(mtx_path, bc_path, ft_path, map_path)) {
    if (!file.exists(f)) stop("Missing input: ", f)
  }
  mat <- readMM(mtx_path)                              # genes x barcodes
  bcs <- readLines(bc_path)                            # barcode per column
  fts <- read.table(ft_path, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE)[[1]]     # gene_id per row
  map <- read.table(map_path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)          # sample_id, barcode

  rownames(mat) <- fts
  colnames(mat) <- bcs

  # keep only barcodes we know are ours (drop any spurious matches)
  keep_bc <- intersect(bcs, map$barcode)
  mat <- mat[, keep_bc, drop = FALSE]

  # rename columns barcode -> sample_id
  bc_to_sid <- setNames(map$sample_id, map$barcode)
  colnames(mat) <- unname(bc_to_sid[colnames(mat)])
  cat(sprintf("Pool %d [%s]: %d genes x %d samples\n",
              p, matrix_name, nrow(mat), ncol(mat)))
  mat
}

merge_pools <- function(matrix_name, out_path) {
  cat("\n=== Merging pools for ", matrix_name, " ===\n", sep = "")
  mats <- lapply(1:4, read_pool, matrix_name = matrix_name)

  # sanity: same gene set across pools?
  gene_ok <- all(vapply(mats[-1], function(m)
                        identical(rownames(m), rownames(mats[[1]])), logical(1)))
  if (!gene_ok) stop("Gene ID sets differ across pools -- can't merge")

  combined <- do.call(cbind, mats)
  cat(sprintf("Combined: %d genes x %d samples\n",
              nrow(combined), ncol(combined)))

  # write TSV in Zea_mays_counts.txt shape (rownames = gene_id, header = sample_ids)
  df <- as.data.frame(as.matrix(combined))
  write.table(df, out_path, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = NA)
  cat("Wrote ", out_path, " (", format(file.size(out_path) / 1e6, digits = 3),
      " MB)\n", sep = "")
  invisible(combined)
}

umi <- merge_pools("umiDedup-1MM_Directional.mtx",
                   file.path(out_dir, "Zea_mays_counts.txt"))
raw <- merge_pools("umiDedup-NoDedup.mtx",
                   file.path(out_dir, "Zea_mays_counts_raw.txt"))

# --- one-line library-size comparison so we can spot dedup rate immediately
umi_lib <- colSums(umi); raw_lib <- colSums(raw)
common  <- intersect(names(umi_lib), names(raw_lib))
dedup_rate <- 1 - (umi_lib[common] / raw_lib[common])
cat(sprintf(
  "\nPer-sample PCR-duplicate rate (median %.1f%%, IQR [%.1f%%, %.1f%%])\n",
  100 * median(dedup_rate, na.rm = TRUE),
  100 * quantile(dedup_rate, 0.25, na.rm = TRUE),
  100 * quantile(dedup_rate, 0.75, na.rm = TRUE)))

cat("\nDone.\n")
