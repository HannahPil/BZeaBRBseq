#!/bin/bash
#BSUB -J 06_Counts
#BSUB -q sara
#BSUB -n 1
#BSUB -o logs/06_Counts.%J.out
#BSUB -e logs/06_Counts.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

# activate conda env so Rscript resolves to env's R (with Rsubread installed)
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

Rscript ../scripts/06_featureCounts_Zm.R Zea_mays
