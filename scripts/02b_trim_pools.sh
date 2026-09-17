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
# where POOL_N is 1..4.
#
# Inputs:
#   $poolDir/BZeaBRB{N}_S{N}_L004_R1_001.fastq.gz
#   $poolDir/BZeaBRB{N}_S{N}_L004_R2_001.fastq.gz
#
# Outputs:
#   $baseDir/trimmed/pool_N_R1.fastq.gz     paired R1 (essentially unchanged)
#   $baseDir/trimmed/pool_N_R2.fastq.gz     paired R2 (trimmed)
#   $baseDir/trimmed/pool_N_trim.log        Trimmomatic summary
# ==============================================================================

set -euo pipefail

pool="${1:-}"
if [[ ! "$pool" =~ ^[1-4]$ ]]; then
  echo "usage: $0 <POOL_N>   (POOL_N = 1, 2, 3, or 4)" >&2
  exit 1
fi

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
poolDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez"

R1_in="${poolDir}/BZeaBRB${pool}_S${pool}_L004_R1_001.fastq.gz"
R2_in="${poolDir}/BZeaBRB${pool}_S${pool}_L004_R2_001.fastq.gz"

outDir="${baseDir}/trimmed"
mkdir -p "$outDir"
R1_out="${outDir}/pool_${pool}_R1.fastq.gz"
R2_out="${outDir}/pool_${pool}_R2.fastq.gz"
R1_unp="${outDir}/pool_${pool}_R1.unpaired.fastq.gz"
R2_unp="${outDir}/pool_${pool}_R2.unpaired.fastq.gz"
log="${outDir}/pool_${pool}_trim.log"

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# Trimmomatic ships adapter FASTAs with the install; locate them for ILLUMINACLIP.
# Prefer NexteraPE-PE.fa (BRB-seq uses Nextera adapters via tagmentation).
adapter_dir=$(dirname "$(readlink -f "$(which trimmomatic)")")/../share/trimmomatic*/adapters
adapter_dir=$(ls -d $adapter_dir 2>/dev/null | head -1)
if [ -z "$adapter_dir" ] || [ ! -d "$adapter_dir" ]; then
  # fallback for other conda layouts
  adapter_dir=$(dirname "$(find /usr/local/usrapps/maize/hdpil/hdpil -name 'NexteraPE-PE.fa' 2>/dev/null | head -1)")
fi
adapter_file="${adapter_dir}/NexteraPE-PE.fa"
if [ ! -f "$adapter_file" ]; then
  echo "ERROR: NexteraPE-PE.fa not found under $adapter_dir" >&2
  echo "  Search manually: find /usr/local/usrapps/maize/hdpil/hdpil -name '*.fa' | grep -i adapter" >&2
  exit 1
fi

echo "=== Trimmomatic PE for pool $pool ==="
echo "adapter file:    $adapter_file"
echo "R1 in:           $R1_in  ($(du -h "$R1_in" | cut -f1))"
echo "R2 in:           $R2_in  ($(du -h "$R2_in" | cut -f1))"
echo "R1 out:          $R1_out"
echo "R2 out:          $R2_out"

trimmomatic PE -threads 16 -phred33 \
  "$R1_in" "$R2_in" \
  "$R1_out" "$R1_unp" \
  "$R2_out" "$R2_unp" \
  ILLUMINACLIP:${adapter_file}:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:28 \
  2>&1 | tee "$log"

# clean up the singleton "unpaired" files (STARsolo only needs the paired set)
rm -f "$R1_unp" "$R2_unp"

echo ""
echo "=== pool $pool trim done ==="
ls -la "$R1_out" "$R2_out"
echo ""
echo "Summary lines from log:"
grep -E "Input Read Pairs|Both Surviving" "$log" || tail -5 "$log"
