#!/usr/bin/env Rscript

# ==============================================================================
# 02 -- Prepare STARsolo barcode files (LOCAL or HPC, one-time per sample set)
#
# STARsolo needs two things per pool:
#   1. A cell-barcode whitelist (one 14 nt barcode per line). Every plate uses
#      the same 96 barcodes (one V5B 96-well barcode plate), so one file suffices.
#   2. A barcode -> sample_id map, per pool (pool = plate). Different sample_ids
#      in each pool map to the same 96 physical barcodes, so this file is per-pool.
#
# Writes a map for EVERY plate in data/metadata.csv that has plate_pos filled,
# so new plates need no edit here: add their rows to metadata.csv and rerun.
#
# Inputs:
#   data/metadata.csv                                sample_id, plate_pos, plate
#   data/barcodes.txt                                Name (plate_pos, e.g. A01) <tab> B1 (barcode)
#     Override the path via BARCODES_TXT env var if needed.
#
# Outputs (all under data/starsolo/):
#   barcode_whitelist.txt              96 barcodes, one per line
#   pool_N_barcode_map.tsv             sample_id \t barcode  (one per plate)
#
# Base R only (no tidyverse) so it runs in a minimal env.
# ==============================================================================

data_dir <- "data"
out_dir  <- file.path(data_dir, "starsolo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bc_path <- Sys.getenv("BARCODES_TXT",
                      unset = file.path(data_dir, "barcodes.txt"))
if (!file.exists(bc_path)) {
  stop("barcodes.txt not found at: ", bc_path,
       "\n  Set BARCODES_TXT=<path> or copy the file into the repo.")
}

# --- read inputs ----------------------------------------------------------
metadata <- read.csv(file.path(data_dir, "metadata.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
barcodes <- read.table(bc_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
stopifnot(all(c("Name", "B1") %in% colnames(barcodes)))
colnames(barcodes)[colnames(barcodes) == "Name"] <- "plate_pos"
colnames(barcodes)[colnames(barcodes) == "B1"]   <- "barcode"

cat("Metadata rows:", nrow(metadata), "\n")
cat("Barcode rows: ", nrow(barcodes), "\n")

# --- validate barcodes ---------------------------------------------------
bad_len <- barcodes$barcode[nchar(barcodes$barcode) != 14]
if (length(bad_len) > 0) stop("Non-14 nt barcode(s): ", paste(bad_len, collapse = ", "))
bad_chr <- barcodes$barcode[grepl("[^ACGT]", barcodes$barcode)]
if (length(bad_chr) > 0) stop("Non-ACGT barcode(s): ", paste(bad_chr, collapse = ", "))
dup_bc <- barcodes$barcode[duplicated(barcodes$barcode)]
if (length(dup_bc) > 0) stop("Duplicate barcode(s): ", paste(dup_bc, collapse = ", "))
cat("Barcodes look clean:", nrow(barcodes), "x 14 nt A/C/G/T, no dupes\n")

# --- whitelist (single file, all 96 barcodes) ----------------------------
whitelist_path <- file.path(out_dir, "barcode_whitelist.txt")
writeLines(barcodes$barcode, whitelist_path)
cat("Wrote ", whitelist_path, " (", nrow(barcodes), " barcodes)\n", sep = "")

# --- per-pool barcode -> sample_id maps ----------------------------------
stopifnot(all(c("sample_id", "plate_pos", "plate") %in% colnames(metadata)))
meta <- metadata[!is.na(metadata$plate) &
                 !is.na(metadata$plate_pos) & metadata$plate_pos != "",
                 c("sample_id", "plate_pos", "plate")]
meta$plate_pos <- trimws(meta$plate_pos)
dup_sid <- meta$sample_id[duplicated(meta$sample_id)]
if (length(dup_sid) > 0) stop("Duplicate sample_id(s) in metadata: ",
                              paste(utils::head(dup_sid, 5), collapse = ", "))
mapped <- merge(meta, barcodes, by = "plate_pos", all.x = FALSE)

missing <- setdiff(meta$sample_id, mapped$sample_id)
if (length(missing) > 0) {
  stop(length(missing), " samples in metadata have a plate_pos with no barcode in ",
       bc_path, ": ", paste(utils::head(missing, 5), collapse = ", "))
}

plates <- sort(unique(as.integer(meta$plate)))
cat("Plates with samples:", paste(plates, collapse = " "), "\n")
for (p in plates) {
  pool_df <- mapped[mapped$plate == p, c("sample_id", "barcode")]
  pool_df <- pool_df[order(pool_df$sample_id), , drop = FALSE]
  path <- file.path(out_dir, sprintf("pool_%d_barcode_map.tsv", p))
  write.table(pool_df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Wrote %s (%d samples)\n", path, nrow(pool_df)))
  dup <- pool_df$barcode[duplicated(pool_df$barcode)]
  if (length(dup) > 0) {
    stop("Pool ", p, " has two samples on the same well/barcode: ",
         paste(dup, collapse = ", "))
  }
}

cat("\nDone.\n")
