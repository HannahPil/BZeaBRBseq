#!/bin/bash

# ==============================================================================
# PPL2 coverage (HPC) — per-base coverage across ppl2 / PnsL1 (Zm00001eb012750)
#
# Question this answers: B73-background samples count ~1 read at ppl2 while
# teosinte carriers count ~98 (55-fold, p ~1e-15). Two explanations survive the
# in-silico checks in BZea_PPL2_2026/scripts/11:
#
#   (a) B73 genuinely does not express ppl2
#   (b) B73 DOES express it, but the gene model's 3' end is wrong, so the reads
#       fall outside the feature that BRB-seq counts
#
# Per-base depth separates them. Under (a) B73 has no reads anywhere in the
# locus; under (b) B73 has reads over the gene body but none in the counted
# 3' window.
#
#   ppl2              chr1:42030859-42032683   (+ strand, 8 exons, 1825 bp)
#   upstream nbr      Zm00001eb012740 ends 42030276   (583 bp gap)
#   downstream nbr    Zm00001eb012760 starts 42032708 (- strand, 25 bp gap)
#
# NOTE the downstream neighbour: ppl2 (+) and Zm00001eb012760 (-) are
# convergent and only 25 bp apart. There is almost no room for an extended
# 3' UTR before the neighbour begins, so explanation (b) is already tightly
# constrained -- and any unstranded counting at this locus would be ambiguous.
# The window is padded to include both flanks so the plot shows the neighbours.
#
# These are the LEGACY per-sample BAMs (one per sample, pre-dedup). The Sep 2026
# STARsolo run aligns per POOL and demultiplexes by barcode internally, so it
# has no per-sample BAMs. Alignment and reference are unchanged, and UMI
# collapsing cannot move a read to a different position, so the coverage
# PROFILE is the same question either way.
#
# samtools depth called ONCE with all BAMs -> single tab-sep matrix
# (chr, pos, depth_sample1, depth_sample2, ...). No per-sample intermediates.
#
# Writes into the cloned repo (hannah/BZeaBRBseq/data/) so that
# `git add data/PPL2_depth_matrix.tsv && git commit && git push` from HPC
# ships the output to local via git (no WinSCP).
#
# Output:
#   $repoDir/data/PPL2_depth_matrix.tsv   (3801 rows) x (2 + N samples) columns
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
alignDir="$baseDir/alignments"
outDir="$repoDir/data"
mkdir -p "$outDir"

REGION="chr1:42029800-42033600"

mapfile -t BAMS < <(ls "$alignDir"/*_Aligned.sortedByCoord.out.bam | sort)
n=${#BAMS[@]}
if [ "$n" -eq 0 ]; then
  echo "ERROR: no BAMs found in $alignDir" >&2
  exit 1
fi
echo "Depth over $REGION for $n BAMs..."

ulimit -n 4096 || true

# ---- index any BAMs missing .bai (samtools depth -r needs an index) ----
missing_idx=()
for b in "${BAMS[@]}"; do
  if [ ! -f "${b}.bai" ] && [ ! -f "${b%.bam}.bai" ]; then
    missing_idx+=("$b")
  fi
done
if [ "${#missing_idx[@]}" -gt 0 ]; then
  echo "Indexing ${#missing_idx[@]} BAMs (missing .bai)..."
  printf '%s\n' "${missing_idx[@]}" \
    | xargs -n 1 -P 4 -I{} samtools index -@ 1 "{}"
  echo "Indexing done."
else
  echo "All BAMs already indexed."
fi

out="$outDir/PPL2_depth_matrix.tsv"
{
  printf "chr\tpos"
  for b in "${BAMS[@]}"; do
    s=$(basename "$b" _Aligned.sortedByCoord.out.bam)
    printf "\t%s" "$s"
  done
  printf "\n"
  samtools depth -a -r "$REGION" "${BAMS[@]}"
} > "$out"

rows=$(wc -l < "$out")
cols=$(head -n 1 "$out" | awk -F'\t' '{print NF}')
echo "Done. Wrote $out ($rows rows x $cols cols)"
