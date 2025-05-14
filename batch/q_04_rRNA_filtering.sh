#!/bin/bash
#BSUB -J 04_rRNA
#BSUB -q sara
#BSUB -n 10
#BSUB -o 04_rRNA.%J.out
#BSUB -e 04_rRNA.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

../scripts/04_rRNA_filtering.sh
