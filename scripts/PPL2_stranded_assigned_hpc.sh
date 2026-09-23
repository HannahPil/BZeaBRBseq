#!/bin/bash

# ==============================================================================
# PPL2 coverage by strand AND assignment (HPC)
#
# The plain coverage figure blends everything at the ppl2 / Zm00001eb012760
# junction into one peak, because `samtools depth` ignores both strand and gene
# assignment. This produces the same profile broken down three ways:
#
#   strand  plus / minus                       (SAM FLAG 0x10)
#   class   ppl2 / nbr / unassigned            (STARsolo's own GX tag)
#   group   Teo / B73 at the ppl2 locus        (CB tag -> sample -> genotype)
#
# so the figure can show directly that the plus-strand pile past ppl2's 3' end
# is assigned to nothing, sits exactly where ppl2's assigned signal stops, and
# is present in both genotypes with a Teo excess on top.
#
# Coverage is accumulated by walking each read's CIGAR over the reference, so no
# temporary per-category BAMs are needed -- one pass per pool.
#
# Needs data/PPL2_sample_groups.tsv (pool, barcode, sample_id, group), written
# locally by make_groups.R and committed.
#
# Output:
#   $repoDir/data/PPL2_depth_by_class.tsv
#     pos, strand, class, group, depth   (summed over all samples in that group)
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
soloDir="$baseDir/starsolo"
groups="$repoDir/data/PPL2_sample_groups.tsv"
out="$repoDir/data/PPL2_depth_by_class.tsv"

[ -s "$groups" ] || { echo "ERROR: missing $groups (run make_groups.R and push)" >&2; exit 1; }

LO=42031500
HI=42033400
REGION="chr1:${LO}-${HI}"
PPL2="Zm00001eb012750"
NBR="Zm00001eb012760"

printf "pos\tstrand\tclass\tgroup\tdepth\n" > "$out"

for pool in 1 2 3 4; do
  bam="$soloDir/pool_${pool}/Aligned.sortedByCoord.out.bam"
  [ -s "$bam" ] || { echo "WARNING: no BAM for pool $pool" >&2; continue; }
  [ -f "${bam}.bai" ] || samtools index -@ 4 "$bam"
  echo "pool $pool ..."

  samtools view -F 0x100 -F 0x800 "$bam" "$REGION" \
  | awk -v pool="$pool" -v grpfile="$groups" -v lo="$LO" -v hi="$HI" \
        -v ppl2="$PPL2" -v nbr="$NBR" '
      BEGIN {
        FS = "\t"; OFS = "\t"
        while ((getline line < grpfile) > 0) {
          split(line, f, "\t")
          if (f[1] == "pool") continue                 # header
          if (f[1] == pool) grp[f[2]] = f[4]           # barcode -> group
        }
      }
      {
        gx = "-"; cb = "-"
        for (i = 12; i <= NF; i++) {
          if ($i ~ /^GX:Z:/) gx = substr($i, 6)
          else if ($i ~ /^CB:Z:/) cb = substr($i, 6)
        }
        if (!(cb in grp)) next                          # sample not in our set
        g = grp[cb]
        s = (int($2 / 16) % 2) ? "minus" : "plus"
        c = (gx == ppl2) ? "ppl2" : (gx == nbr) ? "nbr" : "unassigned"

        # walk the CIGAR over the reference to get the covered span
        pos = $4; cig = $6; n = ""
        for (k = 1; k <= length(cig); k++) {
          ch = substr(cig, k, 1)
          if (ch ~ /[0-9]/) { n = n ch; continue }
          len = n + 0; n = ""
          if (ch == "M" || ch == "D" || ch == "N" || ch == "=" || ch == "X") {
            if (ch != "N" && ch != "D") {
              for (p = pos; p < pos + len; p++)
                if (p >= lo && p <= hi) d[p OFS s OFS c OFS g]++
            }
            pos += len
          }
          # I, S, H, P consume no reference
        }
      }
      END { for (k in d) print k, d[k] }
    ' >> "$out"
done

echo
echo "wrote $out ($(($(wc -l < "$out") - 1)) rows)"
echo
echo "=== total covered bases by strand / class / group ==="
awk -F'\t' 'NR>1 { t[$2 FS $3 FS $4] += $5 } END {
  printf "%-6s %-11s %-4s %14s\n", "strand", "class", "grp", "base-depth"
  for (k in t) { split(k, a, FS)
    printf "%-6s %-11s %-4s %14d\n", a[1], a[2], a[3], t[k] }
}' "$out" | sort -k1,1 -k2,2 -k3,3
