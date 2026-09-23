#!/bin/bash
#BSUB -J PPL2_class
#BSUB -q sara
#BSUB -n 4
#BSUB -o logs/PPL2_class.%J.out
#BSUB -e logs/PPL2_class.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00

# samtools is in the pipeline conda env
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

bash ../scripts/PPL2_stranded_assigned_hpc.sh
