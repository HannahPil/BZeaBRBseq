#!/usr/bin/env Rscript

library(Rsubread)

# Define paths
species <- "Zea_mays"
base_path <- "/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
bam_file <- file.path(base_path, "alignments", "PN4_SID294_Aligned.sortedByCoord.out.bam")
annotation_file <- file.path(base_path, species, paste0(species, ".gtf"))

# Run featureCounts
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

# Print result preview
print(head(counts$counts))
cat("Finished test counting for:", bam_file, "\n")
