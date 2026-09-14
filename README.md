# BZeaBRBseq

RNA-seq analysis pipeline for bulk RNA barcoding sequencing (BRB-seq) of the BZea near-isogenic introgression population. Identifies cis- and trans-eQTLs driven by teosinte introgressions into a B73 maize background.

Preprocessing pipeline (`PIPE_00`–`PIPE_08`) developed by Jonathan Ojeda (Buckler Lab), followed by downstream expression, eQTL, and side analyses.

## Project structure

```
BZeaBRBseq/
├── data/                       # small tables + analysis outputs — tracked
│   ├── external/               # reference genome, teogeno file — gitignored
│   ├── processed/              # big pipeline outputs (counts, edgeR) — gitignored
│   └── FBX_*, EQ_*, ...        # analysis-specific data (prefixed) — tracked
├── output/                     # figures, memos, plots — mostly gitignored
├── scripts/                    # HPC pipeline + local R analyses (prefixed)
├── batch/                      # LSF job submission wrappers (matches script prefixes)
│   └── logs/                   # LSF stdout/stderr — gitignored
└── BZeaBRBseq.Rproj            # RStudio project file
```

## Naming convention

Every script and its per-analysis data / batch wrapper carry a **short prefix** identifying the analysis it belongs to. Single-file analyses just use the prefix as the filename (e.g., `WGCNA.R`, `PCA.R`); multi-file analyses use `PREFIX_specific.ext` (e.g., `eQ_cis.R`, `FBX_analysis4_softclip_hpc.sh`).

| Prefix   | Scope |
|----------|-------|
| `PIPE_`  | Main HPC pipeline stages 00–08 (fastq → count matrix) |
| `EXPR`   | Expression heatmaps of candidate genes |
| `DE`     | Genome-wide differential expression (edgeR) |
| `eQ_`    | eQTL scans (cis, trans, global) |
| `SG`     | Single-gene deep dives (B73 vs teo, co-expression) |
| `WGCNA`  | Weighted gene co-expression network analysis |
| `PCA`    | PCA / batch correction |
| `TPM`    | Counts-to-TPM conversion |
| `EXPORT` | Data tables for collaborators |
| `FBX_`   | fbxl1 mapping-bias investigation (Rubén memo series) |
| `REC_`   | One-off recovery (missing samples, etc.) |

Corresponding data files carry the same prefix (`FBX_depth_matrix.tsv`, `FBX_library_sizes.csv`, etc.) so a `ls data/` groups by analysis at a glance. Shared inputs (metadata, gene lists) stay unprefixed.

## HPC pipeline (`PIPE_00` – `PIPE_08`)

Run on the NCSU sara queue via LSF wrappers in `batch/`. Each script assumes `baseDir=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah` and reads/writes under that path.

| Stage | Script | Description |
|-------|--------|-------------|
| 00 | `PIPE_00_clean_metadata.sh` | Clean and format sample metadata |
| 02 | `PIPE_02_demultiplex.sh` | BRB-seq demultiplexing of pool fastqs into per-sample fastqs |
| 03 | `PIPE_03_trimming_and_QC.sh` | Adapter trimming and quality control |
| 04 | `PIPE_04_rRNA_filtering.sh` | Remove ribosomal RNA reads |
| 05 | `PIPE_05_STAR_alignment.sh` | Align reads to Zea mays genome with STAR |
| 06 | `PIPE_06_featureCounts_Zm.R` | Quantify gene-level read counts (`strandSpecific = 1`) |
| 07 | `PIPE_07_generate_summary_statistics.sh` | Alignment and mapping summaries |
| 08 | `PIPE_08_trimming_stats.sh` | Trimming statistics |

Submit with e.g. `cd batch && bsub < q_PIPE_05_STAR_alignment.sh`.

## Downstream analyses (local R scripts)

All scripts assume the working directory is the project root and read inputs from `data/`.

| Script | Description |
|--------|-------------|
| `PCA.R` | PCA before and after plate batch correction (limma) |
| `EXPR.R` | Expression heatmap of candidate genes across BZea lines |
| `DE.R` | Genome-wide differential expression (edgeR, taxa vs B73) + gene-specific modeling |
| `eQ_cis.R` | Cis-eQTL scan using MatrixEQTL (introgression genotypes as predictors) |
| `eQ_trans.R` | Trans-eQTL scan by source chromosome + Manhattan plots + single-gene trans scans |
| `eQ_global.R` | Genome-wide eQTL architecture contact map, trans hotspot analysis, locus zoom |
| `SG.R` | Single-gene expression plots (B73 vs teosinte, colored by taxa) + co-expression |
| `WGCNA.R` | Weighted gene co-expression network analysis |
| `EXPORT.R` | Export processed data tables for collaborators |
| `TPM.R` | Convert raw counts to TPM using exon-merged gene lengths |

### Key design notes

- **Plates 1–4 only**: all scripts filter to sequenced plates 1–4.
- **Normalization**: the eQTL scripts (`eQ_*`, `SG.R`) use raw library-size log2 CPM; the PCA / DE scripts (`PCA.R`, `EXPR.R`, `DE.R`) use edgeR TMM-normalized log2 CPM. Intentional; see per-script comments.
- **Covariates**: eQTL scripts control for plate; `DE.R` controls for spatial position (Row + Range).
- **Single-gene focus**: `SG.R` has a `FOCUS_GENE` variable at the top for easy gene-by-gene work.

## Side analyses (`FBX_`, `REC_`, …)

Prefixed scripts and data files that don't sit on the main pipeline path. Current examples:

- **FBX** — fbxl1 mapping-bias investigation with Rubén Rellán-Álvarez (memo exchange under `output/FBX_*_memo.pdf`). Sub-analyses use `FBX_analysisN_*`.
- **REC** — one-off recovery pipeline for missing demultiplex samples.

New side analyses go under a new prefix, following the same rule: single-file → `PREFIX.R`; multi-file → `PREFIX_specific.R`.

## Data layout

`data/` is flat with per-file prefixes; two subfolders hold the large or per-machine files:

```
data/
├── external/                       # gitignored; obtain per host
│   ├── Zea_mays.fasta              # from MaizeGDB
│   ├── Zea_mays.gtf                # from MaizeGDB
│   └── results_list_new_name.rds   # teogeno file (Rubén)
├── processed/                      # gitignored; big regenerable outputs
│   ├── Zea_mays_counts.txt         # from PIPE_06
│   ├── edgeR_log2cpm_TMM_filtered.csv  # from DE.R
│   └── edgeR_results_taxa_plus_space.csv
├── metadata.csv, metadata_all.csv, gene_names.csv, …   # shared (unprefixed)
└── FBX_*.csv, FBX_*.tsv, FBX_*.bed, …                   # side-analysis data (prefixed)
```

### Required inputs

| File | Where it lives | Notes |
|------|----------------|-------|
| `data/external/Zea_mays.fasta` | not in git | Reference genome — MaizeGDB B73 v5 |
| `data/external/Zea_mays.gtf` | not in git | Gene annotation — B73 v5 |
| `data/external/results_list_new_name.rds` | not in git | Teosinte introgression segments per BZea genotype |
| `data/processed/Zea_mays_counts.txt` | not in git | Gene-level count matrix (regenerable by rerunning `PIPE_06`) |
| `data/metadata.csv` | tracked | Sample metadata (sample_id, genotype, taxa, plate, Row, Range) |
| `data/Allelic_series_for_expression.csv` | tracked | Which genotypes carry teosinte at each gene |
| `data/candidate_genes.csv` | tracked | Candidate gene list with categories (FT, targ, GWAS, Fst) |
| `data/gene_names.csv` | tracked | Gene ID → gene name mapping |

## Local ↔ HPC workflow (git-first, no WinSCP)

Both the local Windows machine and the HPC clone the same repo. HPC-side scripts write outputs directly into `hannah/BZeaBRBseq/data/` (git-tracked). Small outputs flow local ↔ HPC via `git push` / `git pull` — no WinSCP for anything under ~500 KB. Big files (fasta, gtf, count matrix) stay under `data/external` / `data/processed` and are gitignored.

**Typical loop for a new HPC analysis:**

```bash
# local: write scripts, commit, push
git add scripts/FBX_analysisN_hpc.sh batch/q_FBX_analysisN.sh
git commit -m "..." && git push

# HPC: pull, submit, wait
cd ~/hannah/BZeaBRBseq && git pull
cd batch && bsub < q_FBX_analysisN.sh
bjobs -w

# HPC: after job finishes, commit outputs and push
cd ~/hannah/BZeaBRBseq
git add data/FBX_*                 # only the small pushed files
git commit -m "FBX analysis N: HPC outputs" && git push

# local: pull, run downstream plot
git pull
Rscript scripts/FBX_analysisN_plot.R
```

## Memos

Memos to collaborators (Quarto markdown → HTML → PDF) live at `output/<PREFIX>_*_memo.md`. Render with:

```bash
bash scripts/render_memo.sh output/FBX_reply_memo.md
```

which produces `.html` and `.pdf` next to the source. HTML is Quarto's default Bootstrap theme; PDF is rendered by headless Chrome, matching the style of memos received from collaborators.

Requires: Quarto (`quarto.org`) and Chrome installed at `C:\Program Files\Google\Chrome\Application\chrome.exe` (edit the script if elsewhere).

## R dependencies

`edgeR`, `limma`, `MatrixEQTL`, `Rsubread`, `rtracklayer`, `GenomicRanges`, `tidyverse`, `patchwork`, `ggrepel`, `scales`, `vroom`
