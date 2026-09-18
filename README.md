# BZeaBRBseq

RNA-seq analysis pipeline for bulk RNA barcoding sequencing (BRB-seq) of the BZea near-isogenic introgression population. Identifies cis- and trans-eQTLs driven by teosinte introgressions into a B73 maize background.

Preprocessing pipeline follows Alithea Genomics' July 2026 MERCURIUS™ BRB-seq data-analysis workflow (STARsolo-based, UMI-collapsed counts as default). The original pipeline from Jonathan Ojeda (Buckler Lab) is preserved under `scripts/legacy/` for reference.

## Project structure

```
BZeaBRBseq/
├── data/                       # small tables + analysis outputs — tracked
│   ├── external/               # reference genome, teogeno file — gitignored
│   ├── processed/              # big pipeline outputs (count matrices) — gitignored
│   ├── starsolo/               # STARsolo whitelist + per-pool barcode maps — tracked
│   ├── barcodes.txt            # plate_pos → 14 nt barcode, per V5B kit — tracked
│   └── FBX_*, ...              # analysis-specific data (prefixed) — tracked
├── output/                     # figures, memos, plots — mostly gitignored
├── scripts/                    # pipeline + local R analyses (prefixed)
│   └── legacy/                 # Buckler-lab pipeline (BRBseqTools + STAR + featureCounts)
├── batch/                      # LSF job wrappers
│   ├── legacy/                 # LSF wrappers for scripts/legacy/
│   └── logs/                   # LSF stdout/stderr — gitignored
└── BZeaBRBseq.Rproj            # RStudio project file
```

## Naming convention

Every downstream analysis (not the main pipeline — see below) carries a **short prefix** identifying it. The rule:

- **Commonly-recognized abbreviations** (WGCNA, PCA, TPM, DE) may stand alone as the filename: `WGCNA.R`, `PCA.R`.
- **Less-recognized prefixes** are followed by a descriptive suffix so a new reader knows what's inside: `SG_single_gene_analysis.R`, `EXPR_expressionAnalysis.R`, `EXPORT_for_collaborators.R`.
- **Multi-file analyses** always use the `PREFIX_specific.ext` form: `eQ_cis.R`, `eQ_trans.R`, `FBX_analysis4_softclip_hpc.sh`.

The **main HPC pipeline** keeps its numbered filenames (`00_clean_metadata.sh`, …, `08_trimming_stats.sh`) — the numbering is the identifier.

| Prefix   | Scope |
|----------|-------|
| (00-08 numbered) | Main HPC pipeline stages (fastq → count matrix) |
| `EXPR_`  | Expression heatmaps of candidate genes |
| `DE`     | Genome-wide differential expression (edgeR) |
| `eQ_`    | eQTL scans (cis, trans, global) |
| `SG_`    | Single-gene deep dives (B73 vs teo, co-expression) |
| `WGCNA`  | Weighted gene co-expression network analysis |
| `PCA`    | PCA / batch correction |
| `TPM`    | Counts-to-TPM conversion |
| `EXPORT_`| Data tables for collaborators |
| `FBX_`   | fbxl1 mapping-bias investigation (Rubén memo series) |
| `REC_`   | One-off recovery (missing samples, etc.) |

Corresponding data files carry the same prefix (`FBX_depth_matrix.tsv`, `FBX_library_sizes.csv`, etc.) so a `ls data/` groups by analysis at a glance. Shared inputs (metadata, gene lists) stay unprefixed.

## HPC pipeline (STARsolo, per Alithea July 2026 workflow)

Detailed usage — cold start, adding new samples, rerunning one pool, troubleshooting — lives in [`docs/PIPELINE.md`](docs/PIPELINE.md). Section below is an overview.

Runs on the NCSU sara queue via LSF wrappers in `batch/`. Each script assumes `baseDir=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah`, `repoDir=$baseDir/BZeaBRBseq`, and the STARsolo command comes straight from Alithea's July 2026 data-analysis manual §1.4 — one alignment step handles sample-barcode demux (from R1 first 14 nt), UMI extraction (R1 nt 15–28), adapter clipping, alignment, and UMI-collapsed gene counting. Dual-dedup mode (`--soloUMIdedup "1MM_Directional NoDedup"`) emits both the UMI-collapsed and raw-read matrices in the same run.

| Stage | Script | Description | Where it runs |
|-------|--------|-------------|---------------|
| 02  | `02_prepare_barcodes.R` | Build STARsolo barcode whitelist + per-pool barcode ↔ sample_id maps from `data/metadata.csv` + `data/barcodes.txt`. One-time (until new samples). | LOCAL or HPC |
| 02b | `02b_trim_pools.sh POOL_N` | Trimmomatic PE on pool R1+R2 to trim the 150 PE overshoot before STARsolo (needed because Alithea's STARsolo command is tuned for ~90 nt R2, not 150 nt). Keeps R1/R2 in sync. | HPC (LSF array) |
| 03  | `03_STARsolo_per_pool.sh POOL_N` | STARsolo alignment + demux + UMI count for one pool. Reads trimmed pool fastqs, writes to `hannah/starsolo/pool_N/Solo.out/Gene/raw/`. | HPC (LSF array) |
| 04  | `04_merge_pools.R` | Merge per-pool `.mtx` outputs into two tab-delimited count matrices in `data/processed/`, plus a `starsolo_pool_qc.csv` summary. | HPC or LOCAL |

**Submit commands** (once `02_prepare_barcodes.R` has been run and its outputs are in `data/starsolo/`):

```bash
cd batch
bsub -J "trim[1-4]" < q_02b_trim.sh          # ~1-6 h per pool depending on size
# wait for trim to finish, then:
bsub -J "STARsolo[1-4]" < q_03_STARsolo.sh   # ~1.5-10 h per pool depending on size
# wait for STARsolo to finish, then:
bsub < q_04_merge_pools.sh                   # ~10 min
```

Pool 1 is by far the largest (~40 GB R1); pools 3–4 finish quickly. LSF array `[1-4]` fans all four pools onto separate hosts in parallel.

**Outputs (canonical):**
- `data/processed/Zea_mays_counts.txt` — UMI-collapsed count matrix (this is what downstream analyses read)
- `data/processed/Zea_mays_counts_raw.txt` — raw-read count matrix (kept for PCR-bias comparison)
- `data/processed/starsolo_pool_qc.csv` — per-pool trim + alignment + dedup summary

**External resources the pipeline reads (not in the repo):**
- `/rsstu/users/r/rrellan/sara/ref/STAR_index/` — lab-shared STAR index, built once from Zea mays B73 v5
- `hannah/Zea_mays/Zea_mays.gtf` — gene annotation (also mirrored in `data/external/`)

### Legacy pipeline (Buckler Lab, BRBseqTools + STAR + featureCounts)

Preserved under `scripts/legacy/` (`00_clean_metadata.sh`, `02_demultiplex.sh`, `03_trimming_and_QC.sh`, `04_rRNA_filtering.sh`, `05_STAR_alignment.sh`, `06_featureCounts_Zm.R`, `07`, `08`) and `batch/legacy/`. Kept for reference — the original count matrix (raw read counts, no UMI dedup) was produced by this pipeline. Not intended to be re-run.

## Downstream analyses (local R scripts)

All scripts assume the working directory is the project root and read inputs from `data/`.

| Script | Description |
|--------|-------------|
| `PCA.R` | PCA before and after plate batch correction (limma) |
| `EXPR_expressionAnalysis.R` | Expression heatmap of candidate genes across BZea lines |
| `DE.R` | Genome-wide differential expression (edgeR, taxa vs B73) + gene-specific modeling |
| `eQ_cis.R` | Cis-eQTL scan using MatrixEQTL (introgression genotypes as predictors) |
| `eQ_trans.R` | Trans-eQTL scan by source chromosome + Manhattan plots + single-gene trans scans |
| `eQ_global.R` | Genome-wide eQTL architecture contact map, trans hotspot analysis, locus zoom |
| `SG_single_gene_analysis.R` | Single-gene expression plots (B73 vs teosinte, colored by taxa) + co-expression |
| `WGCNA.R` | Weighted gene co-expression network analysis |
| `EXPORT_for_collaborators.R` | Export processed data tables for collaborators |
| `TPM.R` | Convert raw counts to TPM using exon-merged gene lengths |

### Key design notes

- **Plates 1–4 only**: all scripts filter to sequenced plates 1–4.
- **Normalization**: the eQTL scripts (`eQ_*`, `SG_*`) use raw library-size log2 CPM; the PCA / DE scripts (`PCA.R`, `EXPR_*`, `DE.R`) use edgeR TMM-normalized log2 CPM. Intentional; see per-script comments.
- **Covariates**: eQTL scripts control for plate; `DE.R` controls for spatial position (Row + Range).
- **Single-gene focus**: `SG_single_gene_analysis.R` has a `FOCUS_GENE` variable at the top for easy gene-by-gene work.

## Side analyses (`FBX_`, `REC_`, …)

Prefixed scripts and data files that don't sit on the main pipeline path. Current examples:

- **FBX** — fbxl1 mapping-bias investigation with Rubén Rellán-Álvarez (memo exchange under `output/FBX_*_memo.pdf`). Sub-analyses use `FBX_analysisN_*`.
- **REC** — one-off recovery pipeline for missing demultiplex samples.

New side analyses go under a new prefix, following the same rule: single-file → `PREFIX.R`; multi-file → `PREFIX_specific.R`.

## Data layout

`data/` is flat with per-file prefixes; subfolders hold the large or per-machine files:

```
data/
├── external/                       # gitignored; obtain per host
│   ├── Zea_mays.fasta              # from MaizeGDB
│   ├── Zea_mays.gtf                # from MaizeGDB
│   └── results_list_new_name.rds   # teogeno file (Rubén)
├── processed/                      # gitignored; big regenerable outputs
│   ├── Zea_mays_counts.txt         # STARsolo UMI-collapsed (canonical)
│   ├── Zea_mays_counts_raw.txt     # STARsolo raw reads (comparison)
│   ├── starsolo_pool_qc.csv        # per-pool trim + alignment + dedup metrics
│   ├── edgeR_log2cpm_TMM_filtered.csv  # from DE.R
│   └── edgeR_results_taxa_plus_space.csv
├── starsolo/                       # tracked; STARsolo pipeline inputs
│   ├── barcode_whitelist.txt       # 96 barcodes, from 02_prepare_barcodes.R
│   └── pool_{1,2,3,4}_barcode_map.tsv  # per-pool sample_id ↔ barcode
├── barcodes.txt                    # plate_pos ↔ 14 nt barcode (V5B kit)
├── metadata.csv, metadata_all.csv, gene_names.csv, …   # shared (unprefixed)
└── FBX_*.csv, FBX_*.tsv, FBX_*.bed, …                   # side-analysis data (prefixed)
```

### Required inputs

| File | Where it lives | Notes |
|------|----------------|-------|
| `data/external/Zea_mays.fasta` | not in git | Reference genome — MaizeGDB B73 v5 |
| `data/external/Zea_mays.gtf` | not in git | Gene annotation — B73 v5 |
| `data/external/results_list_new_name.rds` | not in git | Teosinte introgression segments per BZea genotype |
| `data/processed/Zea_mays_counts.txt` | not in git | UMI-collapsed gene count matrix (regenerable by rerunning the STARsolo pipeline: 02b → 03 → 04) |
| `data/starsolo/barcode_whitelist.txt` | tracked | 96 barcodes for STARsolo demux (from `02_prepare_barcodes.R`) |
| `data/starsolo/pool_N_barcode_map.tsv` | tracked | Per-pool barcode → sample_id map (from `02_prepare_barcodes.R`) |
| `data/barcodes.txt` | tracked | Plate-position → 14 nt barcode for the V5B kit |
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
