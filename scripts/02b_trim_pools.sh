#!/bin/bash

# ==============================================================================
# 02b -- Pool-level R2 trimming for STARsolo (HPC)
#
# Alithea's STARsolo command (--clipAdapterType CellRanger4) assumes Read 2
# is ~60-90 cycles. Our sequencing was 150 PE, so R2 has ~60-90 bp of
# adapter/polyA overshoot past the useful cDNA that STARsolo's built-in
# clipping isn't tuned for. This step trims R2 with Trimmomatic PE so pairs
# stay in sync, before STARsolo demux/alignment.
#
# Design decisions:
#   - PE mode so paired reads stay lockstep even when reads get dropped
#   - ILLUMINACLIP targets Nextera+TruSeq adapters (superset covers BRB-seq)
#   - R1 is barcode+UMI (14+14 nt): MINLEN 28 keeps every read with a full
#     barcode+UMI intact. R1 will rarely lose bases — no Illumina adapter
#     inside the first 28 nt.
#   - SLIDINGWINDOW 4:15 for quality trim on R2 3' end
#   - Only writes the paired-output files; unpaired reads dropped
#
# Usage:
#     bash scripts/02b_trim_pools.sh <POOL_N>
# where POOL_N is any pool (= plate) number listed in data/pool_fastqs.tsv.
#
# Inputs:
#   data/pool_fastqs.tsv    pool <tab> R1 <tab> R2 (absolute Hazel paths).
#                           One row per pool; add rows when a new batch lands.
#                           Override the path with POOL_FASTQS=<file>.
#
# Outputs:
#   $baseDir/trimmed/pool_N_R1.fastq.gz     paired R1 (essentially unchanged)
#   $baseDir/trimmed/pool_N_R2.fastq.gz     paired R2 (trimmed)
#   $baseDir/trimmed/pool_N_trim.log        Trimmomatic summary
#
# Outputs are written under .tmp names and renamed only after Trimmomatic
# exits cleanly, so a killed job never leaves files that look finished.
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="${baseDir}/BZeaBRBseq"
manifest="${POOL_FASTQS:-${repoDir}/data/pool_fastqs.tsv}"
threads="${SLURM_CPUS_PER_TASK:-16}"

pool="${1:-}"
if [[ ! "$pool" =~ ^[0-9]+$ ]]; then
  echo "usage: $0 <POOL_N>   (a pool number listed in $manifest)" >&2
  exit 1
fi

# --- look up this pool's FASTQs -------------------------------------------
[ -f "$manifest" ] || { echo "ERROR: pool manifest not found: $manifest" >&2; exit 1; }
row=$(awk -F'\t' -v p="$pool" 'NR > 1 && $1 == p' "$manifest" | tr -d '\r')
n_rows=$(printf '%s' "$row" | grep -c . || true)
if [ "$n_rows" -ne 1 ]; then
  echo "ERROR: pool $pool has $n_rows rows in $manifest (need exactly 1)." >&2
  echo "  Pools listed: $(awk -F'\t' 'NR > 1 {printf "%s ", $1}' "$manifest")" >&2
  exit 1
fi
R1_in=$(printf '%s' "$row" | cut -f2)
R2_in=$(printf '%s' "$row" | cut -f3)
for f in "$R1_in" "$R2_in"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty FASTQ for pool $pool: $f" >&2; exit 1; }
done
[ "$R1_in" != "$R2_in" ] || { echo "ERROR: R1 and R2 are the same file for pool $pool" >&2; exit 1; }

outDir="${baseDir}/trimmed"
mkdir -p "$outDir"
R1_out="${outDir}/pool_${pool}_R1.fastq.gz"
R2_out="${outDir}/pool_${pool}_R2.fastq.gz"
R1_tmp="${outDir}/pool_${pool}_R1.tmp.fastq.gz"
R2_tmp="${outDir}/pool_${pool}_R2.tmp.fastq.gz"
R1_unp="${outDir}/pool_${pool}_R1.unpaired.fastq.gz"
R2_unp="${outDir}/pool_${pool}_R2.unpaired.fastq.gz"
log="${outDir}/pool_${pool}_trim.log"
rm -f "$R1_tmp" "$R2_tmp" "$R1_unp" "$R2_unp"

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil
command -v trimmomatic >/dev/null || { echo "ERROR: trimmomatic not on PATH" >&2; exit 1; }

# Trimmomatic ships adapter FASTAs with the install. Path is fixed by the
# conda env location; confirmed present via `find` pre-flight.
adapter_file="/usr/local/usrapps/maize/hdpil/hdpil/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"
if [ ! -f "$adapter_file" ]; then
  echo "ERROR: adapter file not found at $adapter_file" >&2
  echo "  Search: find /usr/local/usrapps/maize/hdpil/hdpil -name 'NexteraPE-PE.fa'" >&2
  exit 1
fi

echo "=== Trimmomatic PE for pool $pool ==="
echo "manifest:        $manifest"
echo "adapter file:    $adapter_file"
echo "threads:         $threads"
echo "R1 in:           $R1_in  ($(du -h "$R1_in" | cut -f1))"
echo "R2 in:           $R2_in  ($(du -h "$R2_in" | cut -f1))"
echo "R1 out:          $R1_out"
echo "R2 out:          $R2_out"

trimmomatic PE -threads "$threads" -phred33 \
  "$R1_in" "$R2_in" \
  "$R1_tmp" "$R1_unp" \
  "$R2_tmp" "$R2_unp" \
  ILLUMINACLIP:${adapter_file}:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:28 \
  2>&1 | tee "$log"

# Trimmomatic prints this line only when it finishes the whole input.
grep -q "Both Surviving" "$log" || { echo "ERROR: no 'Both Surviving' line in $log; trim incomplete" >&2; exit 1; }
mv -f "$R1_tmp" "$R1_out"
mv -f "$R2_tmp" "$R2_out"

# clean up the singleton "unpaired" files (STARsolo only needs the paired set)
rm -f "$R1_unp" "$R2_unp"

echo ""
echo "=== pool $pool trim done ==="
ls -la "$R1_out" "$R2_out"
echo ""
echo "Summary lines from log:"
grep -E "Input Read Pairs|Both Surviving" "$log" || tail -5 "$log"
