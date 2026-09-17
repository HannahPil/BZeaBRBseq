#!/bin/bash
#BSUB -J trim
#BSUB -q sara
#BSUB -n 16
#BSUB -o logs/trim_pool%I.%J.out
#BSUB -e logs/trim_pool%I.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 8:00

# Submit as job array of 4 (one per pool):
#     bsub -J "trim[1-4]" < q_02b_trim.sh
# Pool 1 takes longest (~40 GB R1); pools 3-4 finish fast.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

bash ../scripts/02b_trim_pools.sh "${LSB_JOBINDEX}"
