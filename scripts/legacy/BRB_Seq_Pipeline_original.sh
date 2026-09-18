#!/bin/bash

# Author: Jingjing Zhai
# Date: 2023-10-05
# Last Modified: 2024-Jan-03
# Version: 0.3
# Usage: bash pipeline.sh index|QC|mapping|stat|featureCounts|readSaturation|gff2bed
# Description: This is a pipeline for RNA-Seq analysis of 3' RNA-Seq data (BRB-Seq)
# Contact: Jingjing Zhai, jz963@cornell.edu, zhaijingjing603@gmail.com

# print the usage if the number of arguments is not equal to 1
if [ $# -ne 1 ]; then
    echo "Usage: bash pipeline.sh index|QC|mapping|stat|featureCounts|readSaturation|gff2bed"
    exit
fi

if [ "${1}" == "index" ];then
    # create a directory for genome index if it does not exist
    if [ ! -d genomeIndex ]; then
        mkdir -p genomeIndex;
    fi
    cat keyFile.txt | while read line
    do
        genomeFA=$(echo ${line} | awk '{print $2}')
        annotation=$(echo ${line} | awk '{print $3}')
        speciesName=$(echo ${line} | awk '{print $1}')
        mkdir -p genomeIndex/${speciesName}
        STAR --runThreadN 32 --runMode genomeGenerate --genomeDir genomeIndex/${speciesName} --genomeFastaFiles ${genomeFA} --sjdbGTFfile ${annotation} --genomeSAindexNbases 13 &
    done
elif [ "${1}" == "QC" ];then
    if [ ! -d 01_trimmed_reads ]; then
        mkdir -p 01_trimmed_reads;
    fi
    outDIR="01_trimmed_reads/"
    readsDIR="demultiplexed/"
    adapter="TruSeq3-SE.fa"
    cut -f3 metadata.txt | sed '1d' | while read line
    do
        java -jar /usr/local/usrapps/maize/hdpil/hdpil/share/trimmomatic-0.39-2/trimmomatic.jar SE -threads 16 -phred33 ${readsDIR}${line}.fastq.gz ${outDIR}${line}_trimmed.fq ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
        pigz ${outDIR}${line}_trimmed.fq
        fastqc ${outDIR}${line}_trimmed.fq.gz -O ${outDIR}
    done
elif [ "${1}" == "mapping" ];then
    if [ ! -d 02_mapping ]; then
	    mkdir -p 02_mapping;
    fi
    outDIR='02_mapping/'
    inDIR="01_trimmed_reads/"
    cat metadata.txt | sed '1d' | while read line
    do
        speciesName=$(echo ${line} | awk '{print $1}')
        genomeIndex="genomeIndex/"${speciesName}
        platePOS=$(echo ${line} | awk '{print $3}')
        # speciesName is EMPTY, skip this line
        if [ "${speciesName}" == "EMPTY" ]; then
            echo "${platePOS} is EMPTY, skip this line"
            continue
        else
            STAR --readFilesCommand zcat \
            --outFileNamePrefix ${outDIR}${platePOS} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --genomeDir ${genomeIndex} \
            --runThreadN 32 \
            --readFilesIn ${inDIR}${platePOS}_trimmed.fq.gz \
            --twopassMode Basic \
            --limitGenomeGenerateRAM 128000000000
        fi
    done
elif [ "${1}" == "stat" ];then
    echo -e "platePOS\tplateID\trawReads\tcleanReads\tmappedReads\tuniqMappedReads" > summary_statistics.txt
    cat metadata.txt | sed '1d' | while read line
    do
        speciesName=$(echo ${line} | awk '{print $1}')
        if [ "${speciesName}" == "EMPTY" ]; then
            echo "${platePOS} is EMPTY, skip this line"
            continue
        else
            platePOS=$(echo ${line} | awk '{print $3}')
            plateID=$(echo ${line} | awk '{print $2}')
            samtools stats 02_mapping/${platePOS}Aligned.sortedByCoord.out.bam > 02_mapping/${platePOS}Aligned.sortedByCoord.out.bam.stats
            rawReads=$(( $(zcat demultiplexed/${platePOS}.fastq.gz | wc -l) / 4 ))
            cleanReads=$(grep 'Number of input reads' 02_mapping/${platePOS}Log.final.out | awk -F '|' '{print $2}')
            uniqMappedReads=$(grep 'Uniquely mapped reads number' 02_mapping/${platePOS}Log.final.out | awk -F '|' '{print $2}')
            mappedReads=$(grep 'reads mapped:' 02_mapping/${platePOS}Aligned.sortedByCoord.out.bam.stats | awk -F ':' '{print $2}')
            echo -e "${platePOS}\t${plateID}\t${rawReads}\t${cleanReads}\t${mappedReads}\t${uniqMappedReads}" >> summary_statistics.txt
        fi
    done
elif [ "${1}" == "featureCounts" ];then
    if [ ! -d 03_featureCounts ]; then
	    mkdir -p 03_featureCounts;
    fi
    cat metadata.txt | sed '1d' | while read line
    do
        speciesName=$(echo ${line} | awk '{print $1}')
        platePOS=$(echo ${line} | awk '{print $3}')
        if [ "${speciesName}" == "EMPTY" ]; then
            echo "${platePOS} is EMPTY, skip this line"
            continue
        else
            echo "featureCounts for ${platePOS}"
            featureCounts \
            --primary \
            -T 16 \
            -t gene \
            -g ID \
            -a genome_annotation/${speciesName}.gff3 \
            -o 03_featureCounts/${speciesName}_${platePOS}_featureCounts.txt \
            02_mapping/${platePOS}Aligned.sortedByCoord.out.bam 2> 03_featureCounts/${speciesName}_${platePOS}_featureCounts.log
        fi
    done
elif [ "${1}" == "readSaturation" ];then
    if [ ! -d 04_readSaturation ]; then
        mkdir -p 04_readSaturation;
    fi
    # define a function for read saturation, the input is the line for metadata
    function readSaturation(){
        line=${1}
        platePOS=$(echo ${line} | awk '{print $3}')
        speciesName=$(echo ${line} | awk '{print $1}')
        RPKM_saturation.py -i 02_mapping/${platePOS}Aligned.sortedByCoord.out.bam -r genome_annotation/${speciesName}.bed -o 04_readSaturation/${speciesName}_${platePOS}
    }
    # export the function and run with parallel
    export -f readSaturation
    cat metadata.txt | sed '1d' | parallel -j 24 readSaturation {}
elif [ "${1}" == "gff2bed" ];then
    # install agat first
    # conda install -c bioconda agat --experimental-solver=libmamba
    # The expression s/\.gff3$/.bed/ is a substitution command that replaces the .gff3 suffix with .bed. The = characters are used to indicate that the expression should be evaluated as a Perl expression.
    ls genome_annotation/*.gff3 | parallel -j 4 'agat_convert_sp_gff2bed.pl -gff {} -o {= s/\.gff3$/.bed/ =}'
elif [ "${1}" == "map2rRNA" ];then 
    if [ ! -d map2rRNA ]; then
        mkdir -p map2rRNA;
    fi
    wget -O Ath_ncRNA.fa.gz https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz
    gunzip Ath_ncRNA.fa.gz
    # build index for rRNA
    mkdir -p bowtie2_index
    bowtie2-build Ath_ncRNA.fa bowtie2_index/rRNA
    cat metadata.txt | sed '1d' | while read line
    do
        platePOS=$(echo ${line} | awk '{print $3}')
        bowtie2 -p 32 -x bowtie2_index/rRNA -U 01_trimmed_reads/${platePOS}_trimmed.fq.gz -S map2rRNA/${platePOS}.sam 2> map2rRNA/${platePOS}.log
    done
    grep 'overall alignment rate' map2rRNA/*log | sed -e 's/ overall alignment rate//g' -e 's/map2rRNA\///g' -e 's/.log:/\t/' > map2rRNA_percentage.txt
else
    echo "See you later!"
fi
