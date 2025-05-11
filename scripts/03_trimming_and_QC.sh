#!/bin/bash

# Set directories - modify paths to match your current structure
outDIR="/workdir/joo29/03_trimmed/"
readsDIR="/workdir/joo29/raw_reads/"

# Create output directory if it doesn't exist
if [ ! -d "$outDIR" ]; then
    mkdir -p "$outDIR"
fi

# Loop through each .fastq.gz file in readsDIR
for file in "${readsDIR}"*.fastq.gz
do
    # Extract the base name of the file (without path and extension)
    baseName=$(basename "$file" .fastq.gz)

    # Define the path for the trimmed file
    trimmedFile="${outDIR}${baseName}_trimmed.fq.gz"

    # Check if the trimmed file already exists
    if [ -f "$trimmedFile" ]; then
        echo "Trimmed file ${trimmedFile} already exists, skipping..."
        continue  # Skip the rest of the loop and move to the next file
    fi

    # Process the file with Trimmomatic
    java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE -threads 60 -phred33 \
        "$file" "${outDIR}${baseName}_trimmed.fq" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Compress the output file
    pigz "${outDIR}${baseName}_trimmed.fq"

    # Run FastQC on the compressed output file
    fastqc "${outDIR}${baseName}_trimmed.fq.gz" -O "${outDIR}"
done

echo "Trimming and quality control completed."

