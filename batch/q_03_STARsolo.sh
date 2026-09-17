#!/bin/bash
#BSUB -J STARsolo
#BSUB -q sara
#BSUB -n 16
#BSUB -o logs/STARsolo_pool%I.%J.out
#BSUB -e logs/STARsolo_pool%I.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 8:00

# One STARsolo run per pool. Submit as a job array of size 4:
#     bsub -J "STARsolo[1-4]" < q_03_STARsolo.sh
# LSB_JOBINDEX (1..4) picks which pool this task processes.
#
# Pool 1 is ~40 GB R1 and takes the longest (~4-6 h). Pools 3-4 finish faster.
# 8 h walltime is comfortable for all four running in parallel.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

bash ../scripts/03_STARsolo_per_pool.sh "${LSB_JOBINDEX}"
