#!/usr/bin/env Rscript

# ==============================================================================
# FBX Analysis 0 — strand-setting sanity check (Rubén memo §6.2)
#
# featureCounts is currently run with strandSpecific=1 in 06_featureCounts_Zm.R.
# Rerun featureCounts with all three strand settings (0/1/2) on a handful of
# BAMs and report the fraction of reads assigned by each. The winning setting
# has a markedly higher assignment rate than the others.
#
# R version — uses Rsubread::featureCounts (already in the env) instead of the
# subread CLI (which is not installed).
# ==============================================================================

suppressPackageStartupMessages(library(Rsubread))

baseDir  <- "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
alignDir <- file.path(baseDir, "alignments")
gtf      <- file.path(baseDir, "Zea_mays", "Zea_mays.gtf")
outDir   <- file.path(baseDir, "FBX_analyses", "analysis0_strand_check")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# --- pick a handful of BAMs (5 is plenty; assignment fractions stabilise fast)
bams <- list.files(alignDir, pattern = "Aligned\\.sortedByCoord\\.out\\.bam$",
                   full.names = TRUE)[1:5]
bams <- bams[!is.na(bams)]
if (length(bams) == 0) stop("No BAMs found in ", alignDir)

cat("Using", length(bams), "BAMs for strand check:\n")
cat(paste0("  ", bams, "\n"), sep = "")
cat("\n")

# --- run featureCounts three ways ------------------------------------------
results <- list()
for (s in c(0L, 1L, 2L)) {
  cat("=== strandSpecific = ", s, " ===\n", sep = "")
  fc <- featureCounts(
    files          = bams,
    annot.ext      = gtf,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType    = "gene_id",
    isPairedEnd     = FALSE,
    strandSpecific  = s,
    primaryOnly     = TRUE,
    nthreads        = 4
  )
  # save per-strand summary + counts
  write.table(fc$stat,
              file.path(outDir, sprintf("strandtest_s%d.summary", s)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  results[[as.character(s)]] <- fc$stat
  cat("\n")
}

# --- compact comparison table ----------------------------------------------
cat("\n=========================================================================\n")
cat("ASSIGNMENT RATES (per-sample %Assigned; higher = correct strand)\n")
cat("=========================================================================\n")

per_sample_pct <- function(stat_df) {
  # stat_df: first column is "Status", remaining columns are per-sample counts
  bam_cols <- setdiff(colnames(stat_df), "Status")
  tot <- colSums(stat_df[, bam_cols, drop = FALSE])
  assigned <- as.numeric(stat_df[stat_df$Status == "Assigned", bam_cols])
  round(100 * assigned / tot, 2)
}

pct_tbl <- data.frame(
  sample = setdiff(colnames(results[["0"]]), "Status"),
  s0_pct = per_sample_pct(results[["0"]]),
  s1_pct = per_sample_pct(results[["1"]]),
  s2_pct = per_sample_pct(results[["2"]]),
  stringsAsFactors = FALSE
)
pct_tbl$sample <- basename(pct_tbl$sample)
print(pct_tbl, row.names = FALSE)

# Overall winner
means <- c(s0 = mean(pct_tbl$s0_pct),
           s1 = mean(pct_tbl$s1_pct),
           s2 = mean(pct_tbl$s2_pct))
cat("\nMean %Assigned across samples:\n")
print(round(means, 2))
cat("\nWinner: strandSpecific = ",
    substr(names(which.max(means)), 2, 2), "\n", sep = "")

# --- also print full summaries (all statuses) -------------------------------
cat("\n=========================================================================\n")
cat("FULL SUMMARIES\n")
cat("=========================================================================\n")
for (s in c(0L, 1L, 2L)) {
  cat("\n--- strandSpecific =", s, "---\n")
  print(results[[as.character(s)]])
}

cat("\nOutputs in: ", outDir, "\n", sep = "")
