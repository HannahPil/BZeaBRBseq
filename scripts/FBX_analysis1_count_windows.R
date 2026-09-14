#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 1 — count reads in the 4 fbxl1 test windows (memo §6.3)
#
# Uses Rsubread::featureCounts on all sample BAMs with strandSpecific = 1
# (confirmed by Analysis 0). The 4 windows come from
# data/FBX_fbxl1_test_windows.bed and BED coordinates are converted to SAF
# (1-based) here.
#
# Reads the BED from the cloned repo (hannah/BZeaBRBseq/data/) and writes
# outputs BACK INTO the repo (same directory) with an FBX_ prefix. That
# eliminates the WinSCP step: after the job finishes, `git add data/FBX_*.csv
# && git commit && git push` from the HPC ships the outputs to local via git.
#
# Output (in hannah/BZeaBRBseq/data/):
#   FBX_window_counts.csv             4 windows x N samples raw counts
#   FBX_library_sizes.csv             total counts per sample (main matrix)
#   FBX_window_counts_normalized.csv  windows x samples, counts / lib_size * 1e6
# ==============================================================================

suppressPackageStartupMessages(library(Rsubread))

baseDir  <- "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir  <- file.path(baseDir, "BZeaBRBseq")
alignDir <- file.path(baseDir, "alignments")
bedFile  <- file.path(repoDir, "data", "FBX_fbxl1_test_windows.bed")
mainCts  <- file.path(baseDir, "Zea_mays", "Zea_mays_counts.txt")
outDir   <- file.path(repoDir, "data")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(bedFile), file.exists(mainCts), dir.exists(alignDir))

# ---- BED -> SAF ------------------------------------------------------------
bed <- read.table(bedFile, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start_bed0", "end_bed0", "name", "score", "strand")
saf <- data.frame(
  GeneID = bed$name,
  Chr    = bed$chr,
  Start  = bed$start_bed0 + 1L,   # BED 0-based half-open -> SAF 1-based closed
  End    = bed$end_bed0,
  Strand = bed$strand,
  stringsAsFactors = FALSE
)
cat("Windows to count:\n")
print(saf)

# ---- gather BAMs ----------------------------------------------------------
bams <- list.files(alignDir, pattern = "Aligned\\.sortedByCoord\\.out\\.bam$",
                   full.names = TRUE)
if (length(bams) == 0) stop("No BAMs in ", alignDir)
cat("\nCounting", length(bams), "BAMs...\n")

# ---- run featureCounts ----------------------------------------------------
fc <- featureCounts(
  files          = bams,
  annot.ext      = saf,
  isGTFAnnotationFile = FALSE,
  isPairedEnd    = FALSE,
  strandSpecific = 1L,             # confirmed by Analysis 0
  primaryOnly    = TRUE,
  nthreads       = 4
)

# clean sample names (strip _Aligned.sortedByCoord.out.bam suffix)
colnames(fc$counts) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "",
                           basename(colnames(fc$counts)))

cat("\nWindow counts (first 8 samples):\n")
print(fc$counts[, 1:min(8, ncol(fc$counts))])

cat("\nAssignment summary (window-level):\n")
print(fc$stat)

# ---- library sizes from main counts matrix -------------------------------
cat("\nReading main counts matrix for library-size normalization...\n")
main <- read.delim(mainCts, check.names = FALSE, row.names = 1)
lib_size <- colSums(main)
cat("  ", length(lib_size), " samples, median library size = ",
    round(median(lib_size)), "\n", sep = "")

# align sample order
common <- intersect(colnames(fc$counts), names(lib_size))
if (length(common) < ncol(fc$counts)) {
  warning("Some BAMs missing from main counts: ",
          paste(setdiff(colnames(fc$counts), common), collapse = ", "))
}
win_counts <- fc$counts[, common, drop = FALSE]
lib_size   <- lib_size[common]

# normalized (CPM-like) window values
win_norm <- sweep(win_counts, 2, lib_size, "/") * 1e6

# ---- save -----------------------------------------------------------------
write.csv(win_counts,
          file.path(outDir, "FBX_window_counts.csv"),
          row.names = TRUE)
write.csv(data.frame(sample_id = names(lib_size), lib_size = lib_size),
          file.path(outDir, "FBX_library_sizes.csv"),
          row.names = FALSE)
write.csv(win_norm,
          file.path(outDir, "FBX_window_counts_normalized.csv"),
          row.names = TRUE)

cat("\nSaved:\n")
cat("  ", file.path(outDir, "FBX_window_counts.csv"), "\n", sep = "")
cat("  ", file.path(outDir, "FBX_library_sizes.csv"), "\n", sep = "")
cat("  ", file.path(outDir, "FBX_window_counts_normalized.csv"), "\n", sep = "")
cat("\nCopy these back to local data/ for the ratio analysis.\n")
