#!/usr/bin/env Rscript

# ==============================================================================
# 02 -- Prepare STARsolo barcode files (LOCAL or HPC, one-time per sample set)
#
# STARsolo needs two things per pool:
#   1. A cell-barcode whitelist (one 14 nt barcode per line). Every plate uses
#      the same 96 barcode sequences, so one whitelist suffices.
#   2. A barcode -> sample_id map, per pool (pool = plate). This is where the
#      barcode SET matters: the sets put the same 96 sequences in different wells.
#      STARsolo itself only uses the whitelist; the map is applied in
#      04_merge_pools.R, so a map fix needs only a rerun of 04, not of 03.
#
# Barcode sets (well -> barcode tables), chosen by metadata$batch:
#   batch 1 (plates 1-4)          data/barcodes.txt         V5B, barcodes listed across rows
#   batch 2 (plates 5, 7, 8, 9)   data/barcodes_v5D_96.txt  Set D: same sequences listed
#                                                           down columns (lab note: batch 2
#                                                           barcodes switched to Set D)
# Only 2 of 96 wells (A01, H12) carry the same barcode in both sets, so using the
# wrong set gives nearly every sample another well's reads.
#
# Plates whose batch has no entry in BARCODE_SETS (e.g. not yet sequenced) get
# NO map, and any old map for them is deleted, so 03/04 stop on them instead of
# silently using the wrong set. Add the batch to BARCODE_SETS once its set is known.
#
# Inputs:
#   data/metadata.csv            sample_id, plate_pos, plate, batch
#   the barcode-set files above  Name (plate_pos, e.g. A01) <tab> B1 (barcode)
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

BARCODE_SETS <- c("1" = file.path(data_dir, "barcodes.txt"),
                  "2" = file.path(data_dir, "barcodes_v5D_96.txt"))

read_set <- function(path) {
  if (!file.exists(path)) stop("barcode set not found: ", path)
  b <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(all(c("Name", "B1") %in% colnames(b)))
  b <- data.frame(plate_pos = trimws(b$Name), barcode = trimws(b$B1),
                  stringsAsFactors = FALSE)
  bad_len <- b$barcode[nchar(b$barcode) != 14]
  if (length(bad_len) > 0) stop(path, ": non-14 nt barcode(s): ", paste(bad_len, collapse = ", "))
  bad_chr <- b$barcode[grepl("[^ACGT]", b$barcode)]
  if (length(bad_chr) > 0) stop(path, ": non-ACGT barcode(s): ", paste(bad_chr, collapse = ", "))
  if (anyDuplicated(b$barcode))   stop(path, ": duplicate barcode(s)")
  if (anyDuplicated(b$plate_pos)) stop(path, ": duplicate well(s)")
  cat(sprintf("Barcode set %s: %d x 14 nt A/C/G/T, no dupes\n", path, nrow(b)))
  b
}
sets <- lapply(BARCODE_SETS, read_set)

# --- whitelist (single file; every set must hold the same 96 sequences) ---
for (k in names(sets)) {
  if (!setequal(sets[[k]]$barcode, sets[[1]]$barcode))
    stop("barcode set for batch ", k, " has different sequences from batch ",
         names(sets)[1], "; one whitelist no longer suffices")
}
whitelist_path <- file.path(out_dir, "barcode_whitelist.txt")
writeLines(sets[[1]]$barcode, whitelist_path)
cat("Wrote ", whitelist_path, " (", length(sets[[1]]$barcode), " barcodes)\n", sep = "")

# --- per-pool barcode -> sample_id maps ----------------------------------
metadata <- read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE,
                     check.names = FALSE, colClasses = "character",
                     na.strings = character(0))
stopifnot(all(c("sample_id", "plate_pos", "plate", "batch") %in% colnames(metadata)))
meta <- metadata[metadata$plate != "" & metadata$plate_pos != "",
                 c("sample_id", "plate_pos", "plate", "batch")]
meta$plate_pos <- trimws(meta$plate_pos)
dup_sid <- meta$sample_id[duplicated(meta$sample_id)]
if (length(dup_sid) > 0) stop("Duplicate sample_id(s) in metadata: ",
                              paste(utils::head(dup_sid, 5), collapse = ", "))

plates <- sort(unique(as.integer(meta$plate)))
cat("Plates with samples:", paste(plates, collapse = " "), "\n")
for (p in plates) {
  pm   <- meta[meta$plate == as.character(p), ]
  path <- file.path(out_dir, sprintf("pool_%d_barcode_map.tsv", p))
  b    <- unique(pm$batch)
  if (length(b) != 1) stop("plate ", p, " has more than one batch in metadata")
  if (!b %in% names(sets)) {
    had <- file.exists(path)
    if (had) file.remove(path)
    cat(sprintf("Plate %d: batch '%s' has no barcode set; no map written%s\n",
                p, b, if (had) " (old map removed)" else ""))
    next
  }
  pool_df <- merge(pm, sets[[b]], by = "plate_pos")
  if (nrow(pool_df) != nrow(pm))
    stop("plate ", p, ": ", nrow(pm) - nrow(pool_df), " well(s) not in ", BARCODE_SETS[[b]])
  if (anyDuplicated(pool_df$barcode)) stop("Pool ", p, " has two samples on the same well/barcode")
  pool_df <- pool_df[order(pool_df$sample_id), c("sample_id", "barcode")]
  write.table(pool_df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Wrote %s (%d samples, batch %s, %s)\n",
              path, nrow(pool_df), b, basename(BARCODE_SETS[[b]])))
}

cat("\nDone.\n")
