#!/bin/bash

# ==============================================================================
# FBX Analysis 0 — strand-setting sanity check (Rubén memo §6.2)
#
# featureCounts is currently run with -s 1 in scripts/06_featureCounts_Zm.R.
# If that setting is wrong for the BRB-seq R2 library, every count in the
# matrix is reading the wrong strand. This script runs featureCounts with
# all three strand settings on a handful of BAMs and reports the fraction
# of reads assigned by each. The winning setting has a markedly higher
# assignment rate than the others.
#
# Interpretation (from memo):
#   * winning -s matches current setting (-s 1)  -> strand is correct, proceed
#   * -s 2 wins                                   -> matrix is on the wrong strand;
#                                                    the rest of the fbxl1 analyses
#                                                    need their strands flipped
#   * -s 0 wins                                   -> library is unstranded (unlikely
#                                                    for BRB-seq); use -s 0 downstream
# ==============================================================================

set -e

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# ---- paths ------------------------------------------------------------------
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
alignDir="${baseDir}/alignments"
gtf="${baseDir}/Zea_mays/Zea_mays.gtf"
outDir="${baseDir}/FBX_analyses/analysis0_strand_check"
mkdir -p "$outDir"

# ---- pick a handful of BAMs ------------------------------------------------
# 5 BAMs is enough to establish which strand setting wins — assignment
# fractions are near-identical across samples for a given -s value.
mapfile -t bams < <(ls "${alignDir}"/*Aligned.sortedByCoord.out.bam 2>/dev/null | head -5)
if [ "${#bams[@]}" -eq 0 ]; then
    echo "ERROR: no BAMs found in $alignDir"
    exit 1
fi
echo "Using ${#bams[@]} BAM(s) for the strand check:"
printf '  %s\n' "${bams[@]}"
echo ""

# ---- run featureCounts three ways ------------------------------------------
for s in 0 1 2; do
    echo "=== featureCounts -s ${s} ==="
    featureCounts \
        -s "$s" \
        -t exon \
        -g gene_id \
        -a "$gtf" \
        -o "${outDir}/strandtest_s${s}.txt" \
        --primary \
        "${bams[@]}"
    echo ""
done

# ---- report ----------------------------------------------------------------
echo ""
echo "=========================================================================="
echo "ASSIGNMENT SUMMARIES (higher = better)"
echo "=========================================================================="
for s in 0 1 2; do
    echo ""
    echo "--- -s ${s} ---"
    cat "${outDir}/strandtest_s${s}.txt.summary"
done

echo ""
echo "=========================================================================="
echo "COMPACT COMPARISON: Assigned reads only"
echo "=========================================================================="
printf "%-10s %s\n" "strand" "assigned_per_sample"
for s in 0 1 2; do
    line=$(grep -m1 "^Assigned" "${outDir}/strandtest_s${s}.txt.summary" | cut -f2-)
    printf "%-10s %s\n" "-s ${s}" "$line"
done

echo ""
echo "Outputs in: ${outDir}"
