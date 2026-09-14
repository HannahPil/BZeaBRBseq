#!/bin/bash
#BSUB -J FBX_a4_softclip
#BSUB -q sara
#BSUB -n 4
#BSUB -o logs/FBX_a4_softclip.%J.out
#BSUB -e logs/FBX_a4_softclip.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 1:00

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

bash ../scripts/FBX_analysis4_softclip_hpc.sh
