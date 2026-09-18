# STARsolo pipeline — quick-start guide

End-to-end guide for running the BRB-seq preprocessing pipeline
(fastq → UMI-collapsed count matrix) on the NCSU sara queue.

Sections:

1. [Prerequisites](#1-prerequisites) — one-time setup, check before you start
2. [Cold start](#2-cold-start-first-run-ever-or-new-user-on-hpc) — from a clean HPC clone
3. [Adding new samples](#3-adding-new-samples-warm-start) — most common ongoing workflow
4. [Rerunning one pool](#4-rerunning-just-one-pool) — after a script edit, or a specific pool failed
5. [Troubleshooting](#5-troubleshooting) — common failure modes we've hit

---

## 1. Prerequisites

Check each once. Only re-check if something changed since last time.

### 1a. Conda env with STAR, Trimmomatic, samtools, R

```bash
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil
which STAR trimmomatic samtools Rscript
STAR --version    # expect 2.7.11b or newer
```

If any of those `which` calls prints "no X in ...", install it into the env:
```bash
conda install -p /usr/local/usrapps/maize/hdpil/hdpil -c bioconda -c conda-forge <missing_tool>
```

### 1b. STAR index present

```bash
ls -la /rsstu/users/r/rrellan/sara/ref/STAR_index/
```
Should show `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, `chrLength.txt`, etc. If empty or missing files, the index needs to be built (~1 h, one-time). Ask Rubén — this is a lab-shared index; someone else may have already rebuilt it.

### 1c. Raw pool fastqs present

```bash
ls -la /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/BZeaBRB{1,2,3,4}_S{1,2,3,4}_L004_R{1,2}_001.fastq.gz
```
Eight files, each multi-GB. If any are missing you need the sequencer output.

### 1d. Metadata matches raw fastqs

`data/metadata.csv` must have `sample_id`, `plate_pos`, `plate` columns and cover every well you sequenced. Check with:
```bash
awk -F',' 'NR>1 {print $10}' data/metadata.csv | sort | uniq -c
```
Expect ~96 samples per plate for plates 1–4.

### 1e. batch/logs directory

```bash
cd ~/hannah/BZeaBRBseq/batch
ls -d logs || mkdir logs
```

---

## 2. Cold start (first run ever, or new user on HPC)

Assumes: prerequisites all pass, but no pipeline outputs exist yet.

### 2a. Build the STARsolo barcode files (LOCAL or HPC)

Runs in seconds. Only needs to happen once per sample set (redo when new samples get added and metadata changes).

```bash
cd ~/hannah/BZeaBRBseq
Rscript scripts/02_prepare_barcodes.R
```

Expected output:
```
Metadata rows: 1382
Barcode rows:  96
Barcodes look clean: 96 x 14 nt A/C/G/T, no dupes
Wrote data/starsolo/barcode_whitelist.txt (96 barcodes)
Wrote data/starsolo/pool_1_barcode_map.tsv (96 samples)
Wrote data/starsolo/pool_2_barcode_map.tsv (96 samples)
Wrote data/starsolo/pool_3_barcode_map.tsv (96 samples)
Wrote data/starsolo/pool_4_barcode_map.tsv (96 samples)
```

Commit + push so the files travel via git:
```bash
git add data/starsolo/
git commit -m "STARsolo barcode files"
git push
```

### 2b. Trim all four pools (HPC, LSF array)

Trimmomatic PE trims the 150 PE overshoot from R2. All four pools in parallel:

```bash
cd batch
bsub -J "trim[1-4]" < q_02b_trim.sh
bjobs -w
```

Expected: 4 jobs in RUN state, one per pool. Pool 1 is the biggest at 40 GB R1 and takes ~5–6 h; pool 4 finishes in ~1 h.

When all four show DONE:
```bash
# quick per-pool trim summary
for p in 1 2 3 4; do
  echo "=== pool $p ==="
  grep "Input Read Pairs" ~/hannah/trimmed/pool_${p}_trim.log
done
```

Expect "Both Surviving" around 78–85%. If any pool drops below 70% survival, investigate before proceeding (bad quality run, adapter mismatch).

### 2c. STARsolo all four pools (HPC, LSF array)

Do this only after all four trim tasks finished cleanly:
```bash
bsub -J "STARsolo[1-4]" < q_03_STARsolo.sh
bjobs -w
```

Pool 1 is ~6–10 h. Pool 4 is ~1.5–2 h. When all four show DONE:

```bash
# quick per-pool alignment check
for p in 1 2 3 4; do
  echo "=== pool $p ==="
  grep -E "Uniquely mapped reads %|Number of input reads" \
       ~/hannah/starsolo/pool_${p}/Log.final.out
done
```

Expect Uniquely-mapped % in the 60–85% range. If any pool is below 50%, something is wrong (bad index, wrong species, low-quality library).

### 2d. Merge pools + build QC summary

```bash
bsub < q_04_merge_pools.sh
bjobs -w
```

~10 minutes. Produces in `data/processed/`:
- `Zea_mays_counts.txt` — UMI-collapsed count matrix (canonical, this is what downstream analyses read)
- `Zea_mays_counts_raw.txt` — raw-read count matrix (comparison, for PCR-bias eyeballing)
- `starsolo_pool_qc.csv` — per-pool trim + alignment + dedup metrics

Commit + push:
```bash
cd ~/hannah/BZeaBRBseq
git add data/processed/starsolo_pool_qc.csv
git commit -m "STARsolo pipeline: pool QC summary"
git push
```

`.txt` count matrices are ~35 MB and gitignored; either leave them per-machine or `git add -f` if you want them tracked.

---

## 3. Adding new samples (warm start)

Say a new sequencing run produces `BZeaBRB5_S5_L004_R1_001.fastq.gz` and `_R2_001.fastq.gz`, and you've added the 96 new sample_id rows to `data/metadata.csv` (plate=5).

### 3a. Regenerate barcode files

```bash
cd ~/hannah/BZeaBRBseq
Rscript scripts/02_prepare_barcodes.R
```
Expect `Wrote data/starsolo/pool_5_barcode_map.tsv (96 samples)` in addition to pools 1–4.

Update the pool count in `02b_trim_pools.sh` and `03_STARsolo_per_pool.sh` — search for `^[1-4]$` in both scripts and change to `^[1-5]$` (or edit to a wider range if you're planning a lot more pools). Commit + push.

### 3b. Trim + align only the new pool

```bash
cd batch
bsub -J "trim[5]" < q_02b_trim.sh          # ~1 h if similar size to pool 4
# after trim[5] finishes:
bsub -J "STARsolo[5]" < q_03_STARsolo.sh   # ~2 h
```

### 3c. Rebuild the merged count matrix

`04_merge_pools.R` needs updating for the new pool. Search the R script for `1:4` and change to `1:5`. Commit + push. Then:
```bash
bsub < q_04_merge_pools.sh
```

Produces the same three files as before, now covering all five pools.

---

## 4. Rerunning just one pool

Common case: a script edit, or one specific pool failed and you want to redo just that one without redoing the others.

```bash
# rerun trim for just pool 2
bsub -J "trim[2]" < q_02b_trim.sh

# rerun STARsolo for just pool 3
bsub -J "STARsolo[3]" < q_03_STARsolo.sh
```

The output directory is per-pool (`hannah/trimmed/pool_2_R*.fastq.gz`, `hannah/starsolo/pool_3/`), so re-running one pool only overwrites that pool's outputs — other pools' outputs stay intact.

Then rerun the merge to fold the updated pool back into the combined count matrix:
```bash
bsub < q_04_merge_pools.sh
```

---

## 5. Troubleshooting

### "Job exited with code 2 in ~13 seconds, empty .err, mostly empty .out"

The script died in its header before printing anything. Almost always a bash `set -e` or `pipefail` failure on:
- A hardcoded path that doesn't exist (`ls -la <path>` to confirm)
- A missing binary (`which <tool>` in the activated conda env)
- A pipe like `ls | head` where the leading command fails on a nonexistent path

Fix: run the script with trace on an **interactive node** (NCSU policy: never on login):
```bash
bsub -Is -q sara -n 2 -W 0:30 bash
# inside interactive shell:
cd ~/hannah/BZeaBRBseq
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil
bash -x scripts/02b_trim_pools.sh 4 2>&1 | head -40
```

The line right before it silently exits is the culprit.

### "EXITING because of FATAL ERROR: could not open genome file"

STAR index is missing or in a different location than the script points to. Check:
```bash
ls -la /rsstu/users/r/rrellan/sara/ref/STAR_index/
```
Should list `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, etc. If the dir is empty or the files are missing, the index was wiped and needs to be rebuilt. See prerequisite 1b.

### "umi_tools dedup returns unpaired garbage" or "STARsolo doesn't produce Solo.out"

The R1 barcode/UMI region isn't where STARsolo expects it. Check that R1 is really the barcode read (not R2):
```bash
zcat /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/BZeaBRB4_S4_L004_R1_001.fastq.gz \
  | head -2
```
Second line should be a 150 nt sequence where the first 28 nt is the barcode+UMI. Should NOT look like a cDNA sequence (would suggest R1/R2 are swapped in the readFilesIn line of `03_STARsolo_per_pool.sh`).

### Trim survival rate <70%

Something's off with adapter or quality. Check the trim log for adapter matches:
```bash
grep "Using" ~/hannah/trimmed/pool_4_trim.log
```
Should list ILLUMINACLIP prefix pairs + Nextera clipping sequences. If the sequences printed don't match Nextera, the wrong adapter file was picked up — update the path in `02b_trim_pools.sh`.

### "git push rejected — remote contains work that you do not have"

HPC pushed something between your local edit and your local push. Just:
```bash
git pull --rebase origin master
git push origin master
```

### Job needs more time than the walltime allowed

Edit the `#BSUB -W` line in the relevant `q_*.sh` wrapper (format is `HH:MM`, e.g. `-W 10:00` for 10 hours), commit, push, `git pull` on HPC, resubmit.

### Killed a job by accident (Ctrl+C during copy-paste in MobaXterm)

Ctrl+C in a terminal SIGINTs the running command. Copy in MobaXterm is `Ctrl+Insert` (or just select text with mouse — auto-copies by default). To restart: same `bsub` command.
