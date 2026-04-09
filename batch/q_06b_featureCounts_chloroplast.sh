#!/bin/bash
#BSUB -J 06b_Counts_cp
#BSUB -q sara
#BSUB -n 1
#BSUB -o 06b_Counts_cp.%J.out
#BSUB -e 06b_Counts_cp.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00

Rscript ../scripts/06b_featureCounts_chloroplast.R
