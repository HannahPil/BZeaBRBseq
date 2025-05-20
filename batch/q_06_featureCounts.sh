#!/bin/bash
#BSUB -J 06_Counts
#BSUB -q sara
#BSUB -n 1
#BSUB -o 06_Counts.%J.out
#BSUB -e 06_Counts.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

Rscript ./06_featureCounts_Zm.R
