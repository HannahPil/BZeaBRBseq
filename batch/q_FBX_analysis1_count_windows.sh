#!/bin/bash
#BSUB -J FBX_a1_count
#BSUB -q sara
#BSUB -n 4
#BSUB -o logs/FBX_a1_count.%J.out
#BSUB -e logs/FBX_a1_count.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 1:00

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

Rscript ../scripts/FBX_analysis1_count_windows.R
