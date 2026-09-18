#!/bin/bash

# ==============================================================================
# PPL2 ortholog check (HPC) — is Zm00001eb012750 really maize PnsL1/PPL2?
#
# Why this matters: in data/NDH_complex_genes_B73v5_COMPLETE.tsv, PnsL1 is the
# weakest call in the table. It rests on "BLAST 63.5% id" to a "PsbP C-terminal
# domain-containing protein", with no GFF3 annotation supporting it. Compare
# CRR3, which got in on 40.3% identity PLUS a GFF3 annotation. PsbP is a large
# family in plants -- PsbP itself, PPL1, PPL2, PPD1-6 -- so a 63.5% hit can
# easily land on the wrong family member.
#
# This matters because PnsL1 is the ONLY single-copy NDH subunit that B73 does
# not express (1.4 log2cpm; the other eleven singletons run 6.7-10.8). If
# Zm00001eb012750 is not the true ortholog, that observation dissolves: it
# would just be a quiet paralog, and the cis-eQTL is real but means something
# different.
#
# Test: reciprocal best hit. Arabidopsis PPL2 (AT2G39470) -> B73v5 proteome,
# then the best maize hit back against the Arabidopsis proteome. PPL1 and PsbP1
# go along as family controls -- if all three land on the same maize gene, the
# family is too collapsed for BLAST alone to resolve and the call needs
# phylogeny instead.
#
# Outputs (written into the repo so git ships them to local):
#   $repoDir/data/PPL2_ortholog_blast_fwd.tsv    AtPPL2/PPL1/PsbP1 -> maize
#   $repoDir/data/PPL2_ortholog_blast_rev.tsv    best maize hits -> Arabidopsis
#   $repoDir/data/PPL2_psbp_family_maize.txt     all maize hits worth checking
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
refDir="$baseDir/reference"
outDir="$repoDir/data"
mkdir -p "$refDir" "$outDir"

MAIZE_PROT="$refDir/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa"
ATH_PROT="$refDir/Athaliana_TAIR10_pep.fa"

# ---- 1. references (fetched beforehand, on the login node) ---------------
# Hazel compute nodes have no outbound internet -- a job that curls here exits
# 28 after ~2 min. scripts/PPL2_ortholog_fetch.sh downloads all three on the
# login node; this script only reads them.
QUERY="$outDir/PPL2_query.faa"
missing=()
[ -s "$MAIZE_PROT" ] || missing+=("$MAIZE_PROT")
[ -s "$ATH_PROT" ]   || missing+=("$ATH_PROT")
[ -s "$QUERY" ]      || missing+=("$QUERY")
if [ "${#missing[@]}" -gt 0 ]; then
  echo "ERROR: missing reference files:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  echo "Run this on the LOGIN node first:" >&2
  echo "  bash scripts/PPL2_ortholog_fetch.sh" >&2
  exit 1
fi
echo "Maize proteome:       $(grep -c '^>' "$MAIZE_PROT") sequences"
echo "Arabidopsis proteome: $(grep -c '^>' "$ATH_PROT") sequences"
echo "Query:                $(grep -c '^>' "$QUERY") sequences"
grep '^>' "$QUERY"

# ---- 3. forward BLAST: Arabidopsis -> maize ------------------------------
FMT="6 qseqid sseqid pident length qcovs evalue bitscore"
if [ ! -s "${MAIZE_PROT}.pdb" ] && [ ! -s "${MAIZE_PROT}.phr" ]; then
  makeblastdb -in "$MAIZE_PROT" -dbtype prot -out "$MAIZE_PROT" >/dev/null
fi
blastp -query "$QUERY" -db "$MAIZE_PROT" -outfmt "$FMT" \
       -max_target_seqs 20 -evalue 1e-5 -num_threads 4 \
       > "$outDir/PPL2_ortholog_blast_fwd.tsv"
echo "Forward hits: $(wc -l < "$outDir/PPL2_ortholog_blast_fwd.tsv")"

echo
echo "=== top maize hits per Arabidopsis query ==="
awk -F'\t' '{if (!seen[$1]++) print $1"\t"$2"\t"$3"% id\tqcov "$5"\tE="$6}' \
    "$outDir/PPL2_ortholog_blast_fwd.tsv"

echo
echo "=== does Zm00001eb012750 appear at all, and where? ==="
grep -n "Zm00001eb012750" "$outDir/PPL2_ortholog_blast_fwd.tsv" || \
  echo "  NOT among the top 20 hits for any query -- the table's call is wrong"

# ---- 4. reciprocal BLAST: best maize hits -> Arabidopsis -----------------
cut -f2 "$outDir/PPL2_ortholog_blast_fwd.tsv" | sort -u > "$outDir/.maize_hits.txt"
# add the annotated gene explicitly so it is tested even if step 3 missed it
grep -o "Zm00001eb012750[^ ]*" "$MAIZE_PROT" | sort -u >> "$outDir/.maize_hits.txt"
sort -u "$outDir/.maize_hits.txt" > "$outDir/PPL2_psbp_family_maize.txt"
rm -f "$outDir/.maize_hits.txt"

seqtk subseq "$MAIZE_PROT" "$outDir/PPL2_psbp_family_maize.txt" \
  > "$outDir/.maize_hits.faa" 2>/dev/null || \
  awk 'NR==FNR{want[$1];next} /^>/{k=substr($1,2); p=(k in want)} p' \
      "$outDir/PPL2_psbp_family_maize.txt" "$MAIZE_PROT" > "$outDir/.maize_hits.faa"

if [ ! -s "${ATH_PROT}.pdb" ] && [ ! -s "${ATH_PROT}.phr" ]; then
  makeblastdb -in "$ATH_PROT" -dbtype prot -out "$ATH_PROT" >/dev/null
fi
blastp -query "$outDir/.maize_hits.faa" -db "$ATH_PROT" -outfmt "$FMT" \
       -max_target_seqs 5 -evalue 1e-5 -num_threads 4 \
       > "$outDir/PPL2_ortholog_blast_rev.tsv"
rm -f "$outDir/.maize_hits.faa"

echo
echo "=== reciprocal best hit for each maize candidate ==="
awk -F'\t' '{if (!seen[$1]++) print $1"\t-> "$2"\t"$3"% id\tE="$6}' \
    "$outDir/PPL2_ortholog_blast_rev.tsv"

echo
echo "Reciprocal best hit holds only if AtPPL2's best maize hit points back to"
echo "AtPPL2. Anything else means the table's PnsL1 assignment needs revisiting."
echo "Done."
