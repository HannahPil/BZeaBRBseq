# ==============================================================================
# EXPORT FOR COLLABORATORS
# Produces aligned CSVs ready to share for ML on the ZeaL population:
#   1. sample_metadata.csv         (sample-level labels for ML targets and covariates)
#   2. expression_counts.csv       (genes x samples, raw featureCounts integers)
#   3. gene_origin_matrix.csv      (genes x samples, 0 = B73, 1 = teosinte introgression)
#   4. introgression_segments.csv  (long format of skim-WGS segments per sample)
#   5. phenotypes_2023.csv         (field phenotypes, spatially corrected; if corr_D4.csv exists)
#   6. phenotypes_2025.csv         (field phenotypes, spatially corrected; if corr_B5.csv exists)
# Outputs: ./output/for_collaborators/
# ==============================================================================

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# ==============================================================================
# 0. OUTPUT DIR
# ==============================================================================

out_dir <- file.path("output", "for_collaborators")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

print("1. Loading data...")

data_dir <- "data"
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
counts   <- read.delim(file.path(data_dir, "Zea_mays_counts.txt"), check.names = FALSE, row.names = 1)
metadata <- read.csv(file.path(data_dir, "inv4m_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)

teogeno <- teogeno[!duplicated(names(teogeno))]

# ==============================================================================
# 2. ALIGN SAMPLES
# ==============================================================================

print("2. Aligning samples...")

sample_df <- metadata |>
  filter(
    sample_id %in% colnames(counts),
    plate %in% c(1, 2, 3, 4)
  ) |>
  mutate(
    genotype_teogeno_key = if_else(genotype == "B73", "B73.B", paste0(genotype, ".B"))
  ) |>
  filter(genotype_teogeno_key %in% names(teogeno)) |>
  arrange(sample_id)

n_in     <- nrow(metadata)
n_kept   <- nrow(sample_df)
n_dropped <- n_in - n_kept
print(paste("   > samples kept:", n_kept, "of", n_in, "  (dropped:", n_dropped, ")"))
print("   > kept per taxa:")
print(table(sample_df$taxa))
dropped_df <- metadata[!(metadata$sample_id %in% sample_df$sample_id), ]
print("   > dropped per taxa:")
print(table(dropped_df$taxa))

# ==============================================================================
# 3. PROCESS GTF COORDINATES
# ==============================================================================

print("3. Processing GTF coordinates...")

gtf <- import(file.path(data_dir, "Zea_mays.gtf"))
gtf_df <- as.data.frame(gtf)

gene_coords <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(
    chr = as.character(seqnames[1]),
    start = min(start),
    end   = max(end),
    .groups = "drop"
  ) |>
  filter(gene_id %in% rownames(counts)) |>
  arrange(gene_id)

genes_gr <- GRanges(
  seqnames = gene_coords$chr,
  ranges   = IRanges(gene_coords$start, gene_coords$end),
  gene_id  = gene_coords$gene_id
)

# ==============================================================================
# 4. BUILD GENE-ORIGIN MATRIX (0 = B73, 1 = teosinte introgression)
# ==============================================================================

print("4. Building gene-origin matrix...")

seg_df <- imap_dfr(teogeno, ~{
  .x |>
    filter(V4 == "Introgression") |>
    transmute(
      key = .y,
      chr = V1,
      start = as.integer(V2),
      end   = as.integer(V3)
    )
})

segs_gr <- GRanges(
  seqnames = seg_df$chr,
  ranges   = IRanges(seg_df$start, seg_df$end),
  key      = seg_df$key
)

hits <- findOverlaps(segs_gr, genes_gr)

origin_mat <- matrix(
  0L,
  nrow = nrow(gene_coords),
  ncol = nrow(sample_df),
  dimnames = list(gene_coords$gene_id, sample_df$sample_id)
)

hit_genes <- genes_gr$gene_id[subjectHits(hits)]
hit_keys  <- segs_gr$key[queryHits(hits)]
key_to_samples <- split(sample_df$sample_id, sample_df$genotype_teogeno_key)

for(i in seq_along(hit_genes)) {
  g_id <- hit_genes[i]
  t_key <- hit_keys[i]
  s_ids <- key_to_samples[[t_key]]
  if(!is.null(s_ids)) {
    origin_mat[g_id, as.character(s_ids)] <- 1L
  }
}

n_genes_w_introgression <- sum(rowSums(origin_mat) > 0)
print(paste("   > genes with introgression in >=1 sample:", n_genes_w_introgression, "of", nrow(origin_mat)))

# ==============================================================================
# 5. WRITE GENE-ORIGIN MATRIX (wide, genes x samples)
# ==============================================================================

print("5. Writing gene_origin_matrix.csv...")

origin_df <- data.frame(
  gene_id = rownames(origin_mat),
  origin_mat,
  check.names = FALSE
)
write.csv(origin_df, file.path(out_dir, "gene_origin_matrix.csv"), row.names = FALSE)

# ==============================================================================
# 6. WRITE INTROGRESSION SEGMENTS (long, one row per sample x segment)
# ==============================================================================

print("6. Writing introgression_segments.csv...")

key_to_sample <- sample_df |>
  select(sample_id, genotype_teogeno_key) |>
  distinct()

segments_out <- seg_df |>
  inner_join(
    key_to_sample,
    by = c("key" = "genotype_teogeno_key"),
    relationship = "many-to-many"
  ) |>
  select(sample_id, chr, start, end) |>
  arrange(sample_id, chr, start)

write.csv(segments_out, file.path(out_dir, "introgression_segments.csv"), row.names = FALSE)

# ==============================================================================
# 7. WRITE SAMPLE METADATA
# ==============================================================================

print("7. Writing sample_metadata.csv...")

sample_metadata <- metadata[metadata$sample_id %in% sample_df$sample_id, ]
sample_metadata <- sample_metadata[order(sample_metadata$sample_id), ]
colnames(sample_metadata)[colnames(sample_metadata) == "teo-species"] <- "teo_species"
sample_metadata <- sample_metadata[, !colnames(sample_metadata) %in% "inv4m_genotype"]

write.csv(sample_metadata, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

# ==============================================================================
# 8. WRITE ALIGNED EXPRESSION COUNTS
# ==============================================================================

print("8. Writing expression_counts.csv...")

counts_aligned <- counts[, sample_df$sample_id]
counts_out <- data.frame(
  gene_id = rownames(counts_aligned),
  counts_aligned,
  check.names = FALSE
)
write.csv(counts_out, file.path(out_dir, "expression_counts.csv"), row.names = FALSE)

# ==============================================================================
# 9. WRITE FIELD PHENOTYPE FILES (year-renamed, joined to sample_id, sorted)
# ==============================================================================

print("9. Writing phenotype files (if present)...")

geno_to_sample <- sample_df[, c("sample_id", "genotype")]

redundant_cols <- c(
  "Plot", "inv4m_introgression",
  "Species", "species", "accession", "species_mex", "origin",
  "notes2"
)

export_phenotypes <- function(in_path, out_name) {
  if (!file.exists(in_path)) {
    print(paste("   > skipped:", in_path, "(not found)"))
    return(invisible(NULL))
  }
  pheno <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)
  pheno <- pheno[, !colnames(pheno) %in% redundant_cols]
  colnames(pheno)[colnames(pheno) == "notes1"] <- "exclusion_note"
  joined <- merge(
    geno_to_sample, pheno,
    by.x = "genotype", by.y = "Genotype",
    all = FALSE
  )
  joined$genotype <- NULL
  joined <- joined[, c("sample_id", setdiff(colnames(joined), "sample_id"))]
  joined <- joined[order(joined$sample_id), ]
  write.csv(joined, file.path(out_dir, out_name), row.names = FALSE)
  print(paste("   > wrote", out_name, "(", nrow(joined), "rows )"))
}

export_phenotypes(file.path(data_dir, "corr_D4.csv"), "phenotypes_2023.csv")
export_phenotypes(file.path(data_dir, "corr_B5.csv"), "phenotypes_2025.csv")

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

print("10. Done. Files written:")
print(paste("   -", file.path(out_dir, "gene_origin_matrix.csv"),
            paste0("  (", nrow(origin_mat), " genes x ", ncol(origin_mat), " samples)")))
print(paste("   -", file.path(out_dir, "introgression_segments.csv"),
            paste0("  (", nrow(segments_out), " segments)")))
print(paste("   -", file.path(out_dir, "sample_metadata.csv"),
            paste0("  (", nrow(sample_metadata), " samples)")))
print(paste("   -", file.path(out_dir, "expression_counts.csv"),
            paste0("  (", nrow(counts_out), " genes x ", ncol(counts_aligned), " samples)")))
