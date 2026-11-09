#!/usr/bin/env Rscript
# ==============================================================================
# META -- Rebuild the plate 5-9 rows of data/metadata.csv from the 2026
# extraction sheet.
#
# The original metadata.csv rows for plates 5-9 were planned layouts. The
# plates actually extracted in 2026 differ:
#   - wells are filled down the columns (A01, B01, ... H01, A02), as in the
#     sheet, not across the rows;
#   - plate boundaries shifted (samples moved between plates 6-9);
#   - three wells became blanks and some wells became checks.
# Source of truth: data/extractions_2026_plates5-9.tsv (copied from the
# BZeaBRB-manifest sheet, tab 2026_extractions).
#
# Sample IDs: each plate keeps its own block of SIDs (plate 5 = 385-480,
# 6 = 481-576, 7 = 577-672, 8 = 673-768, 9 = 769-864), assigned in sheet
# (well) order, so PN<plate>_SID<n> stays unique and plate names do not change.
#
# Conventions for non-line wells follow the existing file:
#   B73          accession B73, taxa B73, check yes
#   Purple Check accession "Purple ", taxa Purple, check no
#   BLANK        accession BLANK, taxa BLANK, check blank (new)
# A new column, extraction_note, carries the sheet's notes (commas -> ";" so the
# unquoted CSV stays one field per column).
# ==============================================================================

meta_path  <- "data/metadata.csv"
sheet_path <- "data/extractions_2026_plates5-9.tsv"
plates     <- 5:9
# Sequencing batch: plates 1-4 are batch 1; plates 5, 7, 8, 9 were sequenced
# together as batch 2 (Novogene, delivered 2026-11-09). Plate 6 is not yet
# sequenced, so its batch stays empty.
batch_of   <- c("5" = "2", "6" = "", "7" = "2", "8" = "2", "9" = "2")

meta  <- read.csv(meta_path, colClasses = "character", check.names = FALSE,
                  na.strings = character(0))   # keep literal "NA" cells as they are
sheet <- read.delim(sheet_path, colClasses = "character", comment.char = "#")
stopifnot(all(c("plate", "well", "plot", "genotype", "note") %in% names(sheet)))
if (!"extraction_note" %in% names(meta)) meta$extraction_note <- ""

# accession -> taxa, learned from every line in the current file
lines <- meta[meta$check == "no" & !meta$taxa %in% c("#N/A", "Purple", "") , ]
acc_taxa <- tapply(lines$taxa, lines$accession, function(x) names(sort(table(x), decreasing = TRUE))[1])

new_rows <- list()
for (p in plates) {
  old <- meta[meta$plate == as.character(p), ]
  sid <- sort(as.integer(sub(".*SID", "", old$sample_id)))
  sh  <- sheet[sheet$plate == as.character(p), ]
  if (nrow(sh) != length(sid)) stop("plate ", p, ": ", nrow(sh), " sheet wells vs ", length(sid), " SIDs")
  if (anyDuplicated(sh$well)) stop("plate ", p, " has a duplicated well in the sheet")

  g   <- sh$genotype
  acc <- ifelse(g == "B73", "B73",
         ifelse(g == "Purple Check", "Purple ",
         ifelse(g == "BLANK", "BLANK", sub("_.*", "", g))))
  tax <- ifelse(g == "B73", "B73",
         ifelse(g == "Purple Check", "Purple",
         ifelse(g == "BLANK", "BLANK", unname(acc_taxa[acc]))))
  if (anyNA(tax)) stop("no taxa known for accession(s): ", paste(unique(acc[is.na(tax)]), collapse = ", "))
  chk <- ifelse(g == "B73", "yes", ifelse(g == "BLANK", "blank", "no"))

  date <- unique(old$date)
  if (length(date) != 1) stop("plate ", p, " has more than one date in metadata")

  new_rows[[length(new_rows) + 1]] <- data.frame(
    sample_id = sprintf("PN%d_SID%d", p, sid), species = "Zea_mays",
    accession = acc, taxa = tax, check = chk, genotype = g, tissue = "leaf",
    date = date, plate_pos = sh$well, plate = as.character(p), batch = unname(batch_of[as.character(p)]),
    plot = ifelse(sh$plot == "NA", "", sh$plot), Row = "", Range = "",
    extraction_note = gsub(",", ";", sh$note, fixed = TRUE), check.names = FALSE, stringsAsFactors = FALSE)
}
new_rows <- do.call(rbind, new_rows)

keep <- meta[!meta$plate %in% as.character(plates), ]
out  <- rbind(keep, new_rows[, names(meta)])
out  <- out[order(as.integer(sub(".*SID", "", out$sample_id))), ]
stopifnot(!anyDuplicated(out$sample_id), nrow(out) == nrow(meta))
out_path <- Sys.getenv("META_OUT", unset = meta_path)
write.csv(out, out_path, row.names = FALSE, quote = FALSE, na = "", eol = "\r\n")
cat("Rewrote plates", paste(plates, collapse = ", "), "->", out_path, "\n")
