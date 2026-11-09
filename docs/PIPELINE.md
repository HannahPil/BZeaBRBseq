# STARsolo pipeline — quick-start guide

End-to-end guide for running the BRB-seq preprocessing pipeline
(fastq → UMI-collapsed count matrix) on a Slurm HPC (NCSU Hazel).

Sections:

1. [Prerequisites](#1-prerequisites) — one-time setup, check before you start
2. [Cold start](#2-cold-start-first-run-ever-or-new-user-on-hpc) — from a clean HPC clone
3. [Adding new samples](#3-adding-new-samples-warm-start) — most common ongoing workflow
4. [Rerunning one pool](#4-rerunning-just-one-pool) — after a script edit, or a specific pool failed
5. [Troubleshooting](#5-troubleshooting) — common failure modes

---

## 1. Prerequisites

Check each once. Only re-check if something changed since last time.

### 1a. Conda env with STAR, Trimmomatic, samtools, R

```bash
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/env
which STAR trimmomatic samtools Rscript
STAR --version    # STARsolo needs 2.7.11b or newer
```

If any of those `which` calls prints "no X in ...", install it into the env:
```bash
conda install -p /path/to/env -c bioconda -c conda-forge <missing_tool>
```

### 1b. STAR index present

```bash
ls -la /path/to/STAR_index/
```
Should show `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, `chrLength.txt`, etc. If empty or the key files are missing, the index needs to be built (~1 h, one-time). Build with `STAR --runMode genomeGenerate` from the reference fasta + gtf.

### 1c. Raw pool fastqs listed in `data/pool_fastqs.tsv`

`data/pool_fastqs.tsv` is the only place the pipeline learns where each pool's raw FASTQs are: one row per pool (pool = plate), tab-separated `pool  R1  R2`, absolute Hazel paths. File names can be anything (Illumina `BZeaBRB1_S1_L004_R1_001.fastq.gz`, Novogene `..._1.fq.gz`, etc.). `02b_trim_pools.sh` stops with a clear message if a pool has no row, two rows, or a missing/empty file.

### 1d. Metadata matches raw fastqs

`data/metadata.csv` must have `sample_id`, `plate_pos`, `plate` columns and cover every well you sequenced. Check that plate counts look right for your kit:
```bash
awk -F',' 'NR>1 {print $10}' data/metadata.csv | sort | uniq -c
```

### 1e. batch/logs directory

```bash
cd <repo>/batch
ls -d logs || mkdir logs
```

---

## 2. Cold start (first run ever, or new user on HPC)

Assumes: prerequisites all pass, but no pipeline outputs exist yet. No script hardcodes the number of pools: the pools run are the ones you list in `--array=`, and the merge takes every finished `starsolo/pool_N/`. The wrappers use absolute paths, so `sbatch` works from any directory.

### 2a. Build the STARsolo barcode files (LOCAL or HPC)

Runs in seconds. Only needs to happen once per sample set (redo when new samples get added and metadata changes).

```bash
cd <repo>
Rscript scripts/02_prepare_barcodes.R
```

Expected output (per pool):
```
Wrote data/starsolo/barcode_whitelist.txt (<N> barcodes)
Wrote data/starsolo/pool_1_barcode_map.tsv (<N> samples)
Wrote data/starsolo/pool_2_barcode_map.tsv (<N> samples)
...
```

Commit + push so the files travel via git:
```bash
git add data/starsolo/
git commit -m "STARsolo barcode files"
git push
```

### 2b. Trim all pools (HPC, Slurm array)

Trimmomatic PE trims the R2 overshoot for long-read sequencing (150 PE and up). All pools in parallel:

```bash
cd batch
sbatch --array=1-4 q_02b_trim.sh
squeue -u $USER
```

Runtime scales with pool size. When all pools have left `squeue`:
```bash
# quick per-pool trim summary
for p in 1 2 3 4; do   # list the pools you ran
  echo "=== pool $p ==="
  grep "Input Read Pairs" /path/to/trimmed/pool_${p}_trim.log
done
```

The `Both Surviving` line reports paired-read survival. Watch for pools that survive substantially less than the others — big spread suggests a library-quality difference worth investigating before proceeding.

### 2c. STARsolo all pools (HPC, Slurm array)

Do this only after all trim tasks finished cleanly:
```bash
sbatch --array=1-4 q_03_STARsolo.sh
squeue -u $USER
```

When all pools have left `squeue`:

```bash
# quick per-pool alignment check
for p in 1 2 3 4; do   # list the pools you ran
  echo "=== pool $p ==="
  grep -E "Uniquely mapped reads %|Number of input reads" \
       /path/to/starsolo/pool_${p}/Log.final.out
done
```

Alithea's July 2026 manual notes that uniquely-mapped % is "typically around 60–85% of MERCURIUS™ DRUG-seq libraries" for mammalian data; expect lower for more repetitive genomes. A big drop from other pools (e.g., one at 20% when the rest are 50%) is the red flag, more than the absolute number.

### 2d. Merge pools + build QC summary

```bash
sbatch q_04_merge_pools.sh
squeue -u $USER
```

~10 minutes. Produces in `data/processed/`:
- `Zea_mays_counts.txt` — UMI-collapsed count matrix (canonical, this is what downstream analyses read)
- `Zea_mays_counts_raw.txt` — raw-read count matrix (comparison, for PCR-bias eyeballing)
- `starsolo_pool_qc.csv` — per-pool trim + alignment + dedup metrics

Commit + push:
```bash
cd <repo>
git add data/processed/starsolo_pool_qc.csv
git commit -m "STARsolo pipeline: pool QC summary"
git push
```

Count-matrix `.txt` files are ~30–40 MB and gitignored; either leave them per-machine or `git add -f` if you want them tracked.

---

## 3. Adding new samples (warm start)

Say a new sequencing run delivers one R1/R2 pair per new plate (e.g. plates 5, 7, 8, 9), and `data/metadata.csv` already has those plates' `sample_id`, `plate_pos` and `plate` rows. No script needs editing.

### 3a. Register the FASTQs and regenerate barcode maps

Add one row per new pool to `data/pool_fastqs.tsv` (absolute Hazel paths to R1 and R2), then:
```bash
cd <repo>
Rscript scripts/02_prepare_barcodes.R
```
It writes a map for every plate in `metadata.csv`. Commit + push, then `git pull` on Hazel.

### 3b. Trim + align only the new pools

```bash
sbatch --array=5,7,8,9 <repo>/batch/q_02b_trim.sh
# after all trim tasks finish:
sbatch --array=5,7,8,9 <repo>/batch/q_03_STARsolo.sh
```
`03` refuses a pool whose trim log lacks `Both Surviving` or whose barcode map is missing, so it cannot run on a half-trimmed pool.

### 3c. Rebuild the merged count matrix

```bash
sbatch <repo>/batch/q_04_merge_pools.sh
```
Merges every `starsolo/pool_N/` and stops, naming the pool, if any is incomplete. To merge a subset: `sbatch --export=ALL,POOLS="1 2 3 4 5" <repo>/batch/q_04_merge_pools.sh`.

---

## 4. Rerunning just one pool

Common case: a script edit, or one specific pool failed and you want to redo just that one without redoing the others.

```bash
# rerun trim for just pool N
sbatch --array=N q_02b_trim.sh

# rerun STARsolo for just pool N
sbatch --array=N q_03_STARsolo.sh
```

The output directory is per-pool (`trimmed/pool_N_R*.fastq.gz`, `starsolo/pool_N/`), so re-running one pool only overwrites that pool's outputs — other pools' outputs stay intact.

Then rerun the merge to fold the updated pool back into the combined count matrix:
```bash
sbatch q_04_merge_pools.sh
```

---

## 5. Troubleshooting

### "Job exited with code 2 in ~13 seconds, empty .err, mostly empty .out"

The script died in its header before printing anything. Almost always a bash `set -e` or `pipefail` failure on:
- A hardcoded path that doesn't exist (`ls -la <path>` to confirm)
- A missing binary (`which <tool>` in the activated conda env)
- A pipe like `ls | head` where the leading command fails on a nonexistent path

Fix: run the script with trace on an **interactive node** (never on login — HPC policy at NCSU and most other sites):
```bash
srun --account=maize_cpu --partition=compute_partners --cpus-per-task=2 --time=00:30:00 --pty bash
# inside interactive shell:
cd <repo>
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/env
bash -x scripts/02b_trim_pools.sh <pool_N> 2>&1 | head -40
```

The line right before it silently exits is the culprit.

### "EXITING because of FATAL ERROR: could not open genome file"

STAR index is missing or in a different location than the script points to. Check:
```bash
ls -la /path/to/STAR_index/
```
Should list `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, etc. If the dir is empty or the files are missing, the index needs rebuilding. See prerequisite 1b.

### STARsolo doesn't produce Solo.out, or the counts look scrambled

The R1 barcode/UMI region isn't where STARsolo expects it. Check that R1 is really the barcode read (not R2):
```bash
zcat /path/to/pool_R1.fastq.gz | head -2
```
Second line should be a full-length read where the first 28 nt is barcode(14) + UMI(14). Should NOT look like a cDNA sequence — that would suggest R1/R2 are swapped in the `readFilesIn` line of `03_STARsolo_per_pool.sh`.

### Trim survival rate looks off (much lower than other pools)

Check the trim log for adapter matches:
```bash
grep "Using" /path/to/trimmed/pool_N_trim.log
```
Should list ILLUMINACLIP prefix pairs + Nextera clipping sequences. If the sequences printed don't match your library's adapters (Nextera for standard BRB-seq V5B kit), the wrong adapter file was picked up — update the path in `02b_trim_pools.sh`.

### "git push rejected — remote contains work that you do not have"

Someone (or you from a different machine) pushed something between your local edit and your local push. Just:
```bash
git pull --rebase origin master
git push origin master
```

### Job needs more time than the walltime allowed

Edit the `#SBATCH --time` line in the relevant `q_*.sh` wrapper (format is `HH:MM:SS`, e.g. `--time=10:00:00` for 10 hours), commit, push, `git pull` on HPC, resubmit.

### Killed a job by accident (Ctrl+C during copy-paste in MobaXterm)

Ctrl+C in a terminal SIGINTs the running command. Copy in MobaXterm is `Ctrl+Insert` (or just select text with mouse — auto-copies by default). To restart: same `sbatch` command.
