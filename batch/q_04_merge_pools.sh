#!/bin/bash
#BSUB -J merge_pools
#BSUB -q sara
#BSUB -n 2
#BSUB -o logs/merge_pools.%J.out
#BSUB -e logs/merge_pools.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 0:30

# Run 04_merge_pools.R on HPC so the .mtx files (many MB each) stay on HPC.
# Only the merged Zea_mays_counts{,_raw}.txt land in the repo, ready to
# git add + push.
#
# Assumes 03_STARsolo produced Solo.out/Gene/raw/ under $baseDir/starsolo/pool_N/.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq

# STARSOLO_ROOT env var points 04_merge_pools.R at the STARsolo outputs on
# HPC rather than the local (unpopulated) data/starsolo/pool_N tree.
STARSOLO_ROOT="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/starsolo" \
  Rscript scripts/04_merge_pools.R
