#!/bin/bash
#BSUB -J UMI_dedup
#BSUB -q sara
#BSUB -n 8
#BSUB -o logs/UMI_dedup.%J.out
#BSUB -e logs/UMI_dedup.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 6:00

# Runs both 06b (dedup) and 06c (featureCounts on dedup BAMs) as one job.
# Dedup ~30-90 min, featureCounts ~1 h, so 6 h walltime is comfortable.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

PARALLEL=8 bash ../scripts/06b_UMI_dedup.sh
Rscript ../scripts/06c_featureCounts_UMI.R Zea_mays
