# BZeaBRBseq

RNA-seq analysis pipeline for bulk RNA barcoding sequencing (BRB-seq) of the BZea near-isogenic introgression population. Identifies cis- and trans-eQTLs driven by teosinte introgressions into a B73 maize background.

Preprocessing pipeline (steps 00-08) developed by Jonathan Ojeda (Buckler Lab), followed by downstream expression and eQTL analysis.

## Project structure

```
BZeaBRBseq/
├── data/               # Input data files (gitignored)
├── output/             # Analysis outputs (gitignored)
├── scripts/            # R analysis scripts (run locally)
├── batch/              # SLURM job submission wrappers
└── BZeaBRBseq.Rproj   # RStudio project file
```

## Pipeline overview

### Preprocessing (HPC / batch scripts)

| Step | Script | Description |
|------|--------|-------------|
| 00 | `00_clean_metadata.sh` | Clean and format sample metadata |
| 02 | `02_demultiplex.sh` | BRB-seq demultiplexing of pool fastqs into per-sample fastqs |
| 03 | `03_trimming_and_QC.sh` | Adapter trimming and quality control |
| 04 | `04_rRNA_filtering.sh` | Remove ribosomal RNA reads |
| 05 | `05_STAR_alignment.sh` | Align reads to Zea mays genome with STAR |
| 06 | `06_featureCounts_Zm.R` | Quantify gene-level read counts |
| 07 | `07_generate_summary_statistics.sh` | Alignment and mapping summaries |
| 08 | `08_trimming_stats.sh` | Trimming statistics |

SLURM wrappers for these steps are in `batch/` (e.g., `q_03_trimming_and_qc.sh`).

### Analysis (local R scripts)

All R scripts assume the working directory is the project root and read input files from `data/`.

| Step | Script | Description |
|------|--------|-------------|
| 10  | `10_PCA.R` | PCA before and after plate batch correction (limma) |
| A1  | `A1_expressionAnalysis.R` | Expression heatmap of candidate genes across BZea lines |
| A2  | `A2_genomeWideAnalysis.R` | Genome-wide differential expression (edgeR, taxa vs B73) and gene-specific modeling |
| A3a | `A3a_cis_eQTL.R` | Cis-eQTL scan using MatrixEQTL (introgression genotypes as predictors) |
| A3b | `A3b_trans_eQTL.R` | Trans-eQTL scan by source chromosome, combined Manhattan plots, single-gene trans scans |
| A3c | `A3c_global_eQTL.R` | Genome-wide eQTL architecture contact map, trans hotspot analysis, locus zoom plots |
| A4  | `A4_single_gene_analysis.R` | Single-gene expression plots (B73 vs teosinte, colored by taxa) and co-expression analysis |
| A5  | `A5_WGCNA.R` | Weighted gene co-expression network analysis |
| A6  | `A6_export_for_collaborators.R` | Export processed data tables for collaborators |
| --  | `counts_to_TPM.R` | Convert raw counts to TPM using exon-merged gene lengths |

### Key design notes

- **Plates 1-4 only**: All scripts filter to sequenced plates 1-4.
- **Normalization**: The eQTL scripts (A3a, A3b, A4) use raw library-size log2 CPM, while the PCA/DE scripts (10, A1, A2) use edgeR TMM-normalized log2 CPM. This is intentional (see comments in each script).
- **Covariates**: eQTL scripts control for plate; genome-wide DE (script A2) controls for spatial position (Row + Range).
- **Single-gene focus**: Script A4 has a `FOCUS_GENE` variable at the top for easy gene-by-gene analysis.

## Required input files (in `data/`)

| File | Description |
|------|-------------|
| `Zea_mays_counts.txt` | Gene-level read count matrix (genes x samples) |
| `metadata.csv` | Sample metadata (sample_id, genotype, taxa, plate, Row, Range) |
| `Zea_mays.gtf` | Zea mays gene annotation (NAM v5) |
| `results_list_new_name.rds` | Teosinte introgression segments per BZea genotype |
| `Allelic_series_for_expression.csv` | Which genotypes carry teosinte at each gene |
| `candidate_genes.csv` | Candidate gene list with categories (FT, targ, GWAS, Fst) |
| `gene_names.csv` | Gene ID to gene name mapping |

## R dependencies

`edgeR`, `limma`, `MatrixEQTL`, `rtracklayer`, `GenomicRanges`, `tidyverse`, `ggrepel`, `scales`, `vroom`
