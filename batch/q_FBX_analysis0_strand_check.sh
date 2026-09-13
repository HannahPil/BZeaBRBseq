#!/bin/bash
#BSUB -J FBX_a0_strand
#BSUB -q sara
#BSUB -n 4
#BSUB -o logs/FBX_a0_strand.%J.out
#BSUB -e logs/FBX_a0_strand.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 1:00

# activate conda env so Rscript resolves to env's R (with Rsubread)
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

Rscript ../scripts/FBX_analysis0_strand_check.R
