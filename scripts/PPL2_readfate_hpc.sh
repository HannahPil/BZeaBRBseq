#!/bin/bash

# ==============================================================================
# PPL2 read fate (HPC) — for every read near the ppl2 / Zm00001eb012760 junction,
# record its STRAND and the gene STARsolo actually assigned it to.
#
# Two questions this answers, which strand-blind `samtools depth` cannot:
#
#   1. Of the reads stacked in the peak just past ppl2's 3' end, what fraction
#      are plus strand (ppl2 orientation, i.e. read-through) versus minus strand
#      (the neighbour's own 3' end)?
#   2. For the plus-strand ones, where do they GO in the counting -- assigned to
#      ppl2, assigned to the neighbour, or discarded as unassigned?
#
# This is possible because the STARsolo run wrote GX/GN gene-assignment tags and
# CB sample barcodes into the BAM (see 03_STARsolo_per_pool.sh --outSAMattributes).
# GX is STARsolo's own answer to "which gene did this read count toward", so no
# inference is required. GX = "-" means the read counted toward nothing.
#
# Layout (chr1):
#   ppl2   Zm00001eb012750  42030859-42032683  PLUS strand
#   gap                     42032684-42032707  24 bp
#   nbr    Zm00001eb012760  42032708-42041803  MINUS strand, 3' end at the LEFT
#
# Windows reported separately, because the interesting one is PEAK:
#   PPL2_3P    42032384-42032683   last 300 bp of ppl2 (what BRB-seq counts)
#   GAP        42032684-42032707   between the two genes
#   PEAK       42032708-42032960   the neighbour's 3' end, where both stack up
#   NBR_BODY   42032961-42041803   the rest of the neighbour
#
# Output (written into the repo so git ships it):
#   $repoDir/data/PPL2_read_fate.tsv
#     pool, window, strand, gene_assigned, n_reads, n_umi
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
soloDir="$baseDir/starsolo"
out="$repoDir/data/PPL2_read_fate.tsv"

REGION="chr1:42032300-42041900"

printf "pool\twindow\tstrand\tgene_assigned\tn_reads\tn_umi\n" > "$out"

for pool in 1 2 3 4; do
  bam="$soloDir/pool_${pool}/Aligned.sortedByCoord.out.bam"
  if [ ! -s "$bam" ]; then
    echo "WARNING: no BAM for pool $pool at $bam" >&2
    continue
  fi
  [ -f "${bam}.bai" ] || samtools index -@ 4 "$bam"
  echo "pool $pool ..."

  # FLAG 16 set   -> read on the minus strand
  # FLAG 16 clear -> read on the plus strand
  # GX tag        -> gene STARsolo assigned the read to ("-" when unassigned)
  # UB tag        -> collapsed UMI ("-" when the read was a duplicate/unusable)
  samtools view -F 0x100 -F 0x800 "$bam" "$REGION" \
  | awk -v pool="$pool" '
      BEGIN { FS = "\t"; OFS = "\t" }
      {
        pos = $4
        # int($2/16)%2 tests FLAG bit 0x10 without gawk-only and()
        strand = (int($2 / 16) % 2) ? "minus" : "plus"
        gx = "-"; ub = "-"
        for (i = 12; i <= NF; i++) {
          if ($i ~ /^GX:Z:/) { gx = substr($i, 6) }
          else if ($i ~ /^UB:Z:/) { ub = substr($i, 6) }
        }
        if      (pos >= 42032384 && pos <= 42032683) w = "PPL2_3P"
        else if (pos >= 42032684 && pos <= 42032707) w = "GAP"
        else if (pos >= 42032708 && pos <= 42032960) w = "PEAK"
        else if (pos >  42032960 && pos <= 42041803) w = "NBR_BODY"
        else                                          w = "OUTSIDE"
        key = pool OFS w OFS strand OFS gx
        n[key]++
        if (ub != "-") u[key]++
      }
      END { for (k in n) print k, n[k], (k in u ? u[k] : 0) }
    ' >> "$out"
done

echo
echo "wrote $out ($(($(wc -l < "$out") - 1)) rows)"
echo
echo "=== totals by window / strand / assigned gene ==="
awk -F'\t' 'NR>1 { r[$2 FS $3 FS $4] += $5 } END {
  printf "%-10s %-6s %-18s %10s\n", "window", "strand", "assigned", "reads"
  for (k in r) { split(k, a, FS); printf "%-10s %-6s %-18s %10d\n", a[1], a[2], a[3], r[k] }
}' "$out" | sort -k1,1 -k2,2 -k4,4nr
