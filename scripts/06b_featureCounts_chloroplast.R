#!/usr/bin/env Rscript

# Feature counting for chloroplast-aligned BAMs
# Counts reads per chloroplast gene using the plastome GTF (converted from GFF3)
# Output: chloroplast_counts.txt in the Zea_mays directory

library(Rsubread)

# Set paths
base_path <- "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
bam_directory <- file.path(base_path, "chloroplast_alignments")
annotation_file <- file.path(base_path, "Zea_mays", "KF241981.1_chloroplast.gtf")

# List all BAM files
bam_files <- list.files(path = bam_directory, pattern = "\\.bam$", full.names = TRUE)

# Stop if no BAM files are found
if (length(bam_files) == 0) {
  stop("No BAM files found in: ", bam_directory)
}

# Initialize a data frame to store count data
count_data <- NULL

# Loop through each BAM file and run featureCounts
for (bam_file in bam_files) {
  counts <- featureCounts(
    files = bam_file,
    annot.ext = annotation_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    isPairedEnd = FALSE,
    primaryOnly = TRUE,
    strandSpecific = 1,
    nthreads = 4
  )

  # Extract the sample ID from the filename
  sample_id <- gsub("_Aligned.*", "", basename(bam_file))

  # Add count data
  if (is.null(count_data)) {
    count_data <- counts$counts
    colnames(count_data) <- sample_id
  } else {
    count_data <- cbind(count_data, counts$counts)
    colnames(count_data)[ncol(count_data)] <- sample_id
  }
}

# Add gene IDs as row names
row.names(count_data) <- counts$annotation$GeneID

# Output file
output_file <- file.path(base_path, "Zea_mays", "chloroplast_counts.txt")
write.table(count_data, file = output_file, row.names = TRUE, sep = "\t", quote = FALSE)

cat("Chloroplast counting completed. Output:", output_file, "\n")
cat("Genes counted:", nrow(count_data), "\n")
cat("Samples counted:", ncol(count_data), "\n")
