#!/usr/bin/env Rscript

# ==============================================================================
# 02 -- Prepare STARsolo barcode files (LOCAL, one-time per sample set)
#
# STARsolo needs two things per pool:
#   1. A cell-barcode whitelist (one 14 nt barcode per line). Same 96 barcodes
#      across all four pools (single V5B 96-well plate), so one file suffices.
#   2. A barcode -> sample_id map, per pool. Different sample_ids in each pool
#      map to the same 96 physical barcodes, so this file is per-pool.
#
# Inputs (LOCAL):
#   data/metadata.csv                                sample_id, plate_pos, plate
#   $HPWORKING/barcodes.txt                          plate_pos -> barcode
#     (default path lives on HPC; for local runs the script also accepts an
#      explicit path via BARCODES_TXT env var)
#
# Outputs (LOCAL, then copied to HPC or pulled via git):
#   data/starsolo/barcode_whitelist.txt              96 barcodes, one per line
#   data/starsolo/pool_N_barcode_map.tsv             sample_id \t barcode
#     for N in 1..4, per pool
#
# Also verifies:
#   - every metadata sample has a plate_pos with a barcode
#   - all barcodes are 14 nt A/C/G/T
#   - no duplicate barcodes within a pool
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

data_dir  <- "data"
out_dir   <- file.path(data_dir, "starsolo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# barcodes.txt lives on HPC by default; allow override for local testing
bc_path <- Sys.getenv(
  "BARCODES_TXT",
  unset = "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hpworking/barcodes.txt"
)
if (!file.exists(bc_path)) {
  stop("barcodes.txt not found at: ", bc_path,
       "\n  Set BARCODES_TXT=<path> or copy the file locally.")
}

# --- read inputs ----------------------------------------------------------
metadata <- read.csv(file.path(data_dir, "metadata.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
barcodes <- read.table(bc_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
# barcodes.txt format: Name (plate_pos, e.g. A01), B1 (barcode)
stopifnot(all(c("Name", "B1") %in% colnames(barcodes)))
colnames(barcodes) <- c("plate_pos", "barcode")

cat("Metadata rows: ", nrow(metadata), "\n")
cat("Barcode rows:  ", nrow(barcodes), "\n")

# --- validate barcodes ---------------------------------------------------
bad_len <- barcodes$barcode[nchar(barcodes$barcode) != 14]
if (length(bad_len) > 0) stop("Non-14 nt barcode(s): ", paste(bad_len, collapse = ", "))
bad_chr <- barcodes$barcode[grepl("[^ACGT]", barcodes$barcode)]
if (length(bad_chr) > 0) stop("Non-ACGT barcode(s): ", paste(bad_chr, collapse = ", "))
dup_bc <- barcodes$barcode[duplicated(barcodes$barcode)]
if (length(dup_bc) > 0) stop("Duplicate barcode(s): ", paste(dup_bc, collapse = ", "))
cat("Barcodes look clean: ", nrow(barcodes), " x 14 nt A/C/G/T, no dupes\n")

# --- whitelist (single file, all 96 barcodes) ----------------------------
whitelist_path <- file.path(out_dir, "barcode_whitelist.txt")
writeLines(barcodes$barcode, whitelist_path)
cat("Wrote ", whitelist_path, " (", nrow(barcodes), " barcodes)\n", sep = "")

# --- per-pool barcode -> sample_id maps ----------------------------------
# Join metadata to barcodes by plate_pos, split by pool.
mapped <- metadata |>
  dplyr::select(sample_id, plate_pos, plate) |>
  dplyr::filter(plate %in% 1:4, !is.na(plate_pos), plate_pos != "") |>
  dplyr::inner_join(barcodes, by = "plate_pos")

missing <- setdiff(metadata$sample_id[metadata$plate %in% 1:4],
                   mapped$sample_id)
if (length(missing) > 0) {
  warning(length(missing), " samples in metadata missing barcodes: ",
          paste(head(missing, 5), collapse = ", "), " ...")
}

for (p in 1:4) {
  pool_df <- mapped |>
    dplyr::filter(plate == p) |>
    dplyr::select(sample_id, barcode) |>
    dplyr::arrange(sample_id)
  path <- file.path(out_dir, sprintf("pool_%d_barcode_map.tsv", p))
  write.table(pool_df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Wrote %s (%d samples)\n", path, nrow(pool_df)))

  # sanity: no duplicate barcodes within a pool
  dup <- pool_df$barcode[duplicated(pool_df$barcode)]
  if (length(dup) > 0) {
    warning("Pool ", p, " has duplicate barcodes: ",
            paste(dup, collapse = ", "))
  }
}

cat("\nDone.\n")
