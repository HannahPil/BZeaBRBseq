#!/usr/bin/env Rscript

# ==============================================================================
# 06c -- featureCounts on UMI-deduplicated BAMs
#
# Parallel of 06_featureCounts_Zm.R, but reads BAMs from alignments_dedup/
# and writes to a separate output (Zea_mays_counts_UMI.txt) so the raw-read
# counts stay intact for comparison.
#
# All parameters mirror 06_featureCounts_Zm.R exactly (strandSpecific=1,
# primaryOnly=TRUE, isPairedEnd=FALSE), so any difference in the output
# matrix vs Zea_mays_counts.txt is attributable to UMI dedup alone.
#
# Usage:
#   Rscript ../scripts/06c_featureCounts_UMI.R Zea_mays
# ==============================================================================

suppressPackageStartupMessages(library(Rsubread))

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
if (is.null(species) || species == "") stop("Usage: 06c_featureCounts_UMI.R <species>")

base_path <- "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
bam_directory <- file.path(base_path, "alignments_dedup")
annotation_file <- file.path(base_path, species, paste0(species, ".gtf"))

bam_files <- list.files(path = bam_directory, pattern = "_dedup\\.bam$",
                        full.names = TRUE)
if (length(bam_files) == 0) {
  stop("No *_dedup.bam files in ", bam_directory,
       " -- run 06b_UMI_dedup.sh first")
}
cat("Counting UMI-deduplicated reads across ", length(bam_files), " BAMs...\n", sep = "")

# One call across all BAMs is much faster than the per-BAM loop in
# 06_featureCounts_Zm.R; behavior is identical.
fc <- featureCounts(
  files               = bam_files,
  annot.ext           = annotation_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  isPairedEnd         = FALSE,
  primaryOnly         = TRUE,
  strandSpecific      = 1,
  nthreads            = 4
)

# strip _dedup suffix from column names so they match Zea_mays_counts.txt
colnames(fc$counts) <- sub("_dedup\\.bam$", "",
                           sub("_Aligned\\.sortedByCoord\\.out", "",
                               basename(colnames(fc$counts))))

output_file <- file.path(base_path, species, paste0(species, "_counts_UMI.txt"))
write.table(fc$counts, file = output_file,
            row.names = TRUE, sep = "\t", quote = FALSE)

cat("Wrote UMI-deduplicated counts to: ", output_file, "\n", sep = "")

# Summary: total counts before / after so you get an immediate dedup-rate readout
lib_totals <- colSums(fc$counts)
cat(sprintf("\nUMI-deduplicated library sizes: median = %.2fM, IQR = [%.2fM, %.2fM]\n",
            median(lib_totals) / 1e6,
            quantile(lib_totals, 0.25) / 1e6,
            quantile(lib_totals, 0.75) / 1e6))
cat("(compare against Zea_mays_counts.txt colSums for the raw-read library sizes)\n")
