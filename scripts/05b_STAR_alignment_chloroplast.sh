#!/bin/bash

# Chloroplast genome alignment using STAR
# Aligns clean reads against the maize plastome (KF241981.1)
# Uses the same clean FASTQ files from the nuclear alignment pipeline

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
metadata="${baseDir}/metadata.txt"
inDIR="${baseDir}/clean_reads/"
outDIR="${baseDir}/chloroplast_alignments/"
genomeFasta="${baseDir}/Zea_mays/KF241981.1_chloroplast.fasta"
gff3File="${baseDir}/Zea_mays/KF241981.1_chloroplast.gff3"
gtfFile="${baseDir}/Zea_mays/KF241981.1_chloroplast.gtf"
genomeDir="${baseDir}/Zea_mays/STAR_index_chloroplast/"
sampleIDs="${baseDir}/sample_ids_chloroplast.txt"

# ---------------------------------------------------------------------------
# Step 1: Convert GFF3 to GTF (requires gffread from cufflinks/stringtie)
# ---------------------------------------------------------------------------
if [ ! -f "$gtfFile" ]; then
    echo "Converting GFF3 to GTF..."
    gffread "$gff3File" -T -o "$gtfFile"
    echo "GTF conversion complete: $gtfFile"
else
    echo "GTF already exists, skipping conversion."
fi

# ---------------------------------------------------------------------------
# Step 2: Build STAR index for the plastome
# ---------------------------------------------------------------------------
mkdir -p "$genomeDir"

if [ ! -f "${genomeDir}/SA" ]; then
    echo "Building STAR index for chloroplast genome..."

    # Get genome length to set genomeSAindexNbases appropriately
    # For a ~140 kb genome, min(14, log2(140000)/2 - 1) ≈ 7
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$genomeDir" \
        --genomeFastaFiles "$genomeFasta" \
        --sjdbGTFfile "$gtfFile" \
        --genomeSAindexNbases 7 \
        --runThreadN 12
    echo "STAR index built."
else
    echo "STAR index already exists, skipping."
fi

# ---------------------------------------------------------------------------
# Step 3: Align clean reads to chloroplast genome
# ---------------------------------------------------------------------------
mkdir -p "$outDIR"

# Extract relevant sample IDs (same as nuclear alignment)
awk -v sp="Zea_mays" -F'\t' '$2 == sp {print $1}' "$metadata" > "$sampleIDs"

while read sample_id; do
    fastqFile="${inDIR}${sample_id}_clean.fastq.gz"

    if [ -f "$fastqFile" ]; then
        if [ ! -f "${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
            echo "Running chloroplast alignment for $sample_id..."
            STAR \
                --outFileNamePrefix "${outDIR}${sample_id}_" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --genomeDir "${genomeDir}" \
                --runThreadN 12 \
                --readFilesCommand zcat \
                --readFilesIn "$fastqFile" \
                --twopassMode Basic \
                --alignIntronMax 1 \
                --alignIntronMin 0
            # Index the BAM
            samtools index "${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam"
        else
            echo "Chloroplast alignment already completed for $sample_id, skipping..."
        fi
    else
        echo "FASTQ file not found for $sample_id"
    fi
done < "$sampleIDs"

echo "All chloroplast alignments completed."
