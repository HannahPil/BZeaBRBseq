#!/bin/bash

# ==============================================================================
# PPL2 ortholog check — STEP 1 of 2, run on the LOGIN NODE.
#
# Hazel compute nodes have no outbound internet (a batch job that curls exits
# 28, curl's timeout code, after ~2 min). So every download happens here, on
# the login node, and PPL2_ortholog_check_hpc.sh then runs BLAST offline.
#
# Downloads nothing it already has, so it is safe to re-run.
#
# Fetches:
#   $refDir/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa   (MaizeGDB)
#   $refDir/Athaliana_TAIR10_pep.fa                           (UniProt reviewed)
#   $outDir/PPL2_query.faa        AtPPL2 + PPL1 + PsbP1 (family controls)
#
# Then submit batch/q_PPL2_ortholog.sh.
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
refDir="$baseDir/reference"
outDir="$repoDir/data"
mkdir -p "$refDir" "$outDir"

MAIZE_PROT="$refDir/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa"
ATH_PROT="$refDir/Athaliana_TAIR10_pep.fa"
QUERY="$outDir/PPL2_query.faa"

# ---- 1. B73v5 proteome ---------------------------------------------------
if [ -s "$MAIZE_PROT" ]; then
  echo "maize proteome already present ($(grep -c '^>' "$MAIZE_PROT") seqs)"
else
  echo "fetching B73v5 proteome from MaizeGDB..."
  curl -fSL --connect-timeout 30 -o "${MAIZE_PROT}.gz" \
    "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.gz"
  gunzip -f "${MAIZE_PROT}.gz"
  echo "  got $(grep -c '^>' "$MAIZE_PROT") sequences"
fi

# ---- 2. Arabidopsis reviewed proteome ------------------------------------
if [ -s "$ATH_PROT" ]; then
  echo "arabidopsis proteome already present ($(grep -c '^>' "$ATH_PROT") seqs)"
else
  echo "fetching Arabidopsis reviewed proteome from UniProt..."
  curl -fSL --connect-timeout 30 -o "$ATH_PROT" \
    "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:3702+AND+reviewed:true&format=fasta"
  echo "  got $(grep -c '^>' "$ATH_PROT") sequences"
fi

# ---- 3. query sequences --------------------------------------------------
# Pulled by gene name rather than hard-coded accession, so the sequences are
# whatever UniProt currently holds rather than something transcribed by hand.
if [ -s "$QUERY" ] && [ "$(grep -c '^>' "$QUERY")" -ge 3 ]; then
  echo "query file already present ($(grep -c '^>' "$QUERY") seqs)"
else
  : > "$QUERY"
  for g in PPL2 PPL1 PSBP1; do
    echo "  fetching At $g ..."
    curl -fSL --connect-timeout 30 \
      "https://rest.uniprot.org/uniprotkb/stream?query=gene:${g}+AND+organism_id:3702+AND+reviewed:true&format=fasta" \
      >> "$QUERY"
  done
  echo "  query file: $(grep -c '^>' "$QUERY") sequences"
fi

echo
echo "=== query sequences ==="
grep '^>' "$QUERY"
echo
echo "All references present. Now submit:  cd batch && bsub < q_PPL2_ortholog.sh"
