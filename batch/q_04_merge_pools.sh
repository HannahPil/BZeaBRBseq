#!/bin/bash
#SBATCH --job-name=merge_pools
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/merge_pools.%j.out
#SBATCH --error=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/merge_pools.%j.err

# Run 04_merge_pools.R on HPC so the .mtx files (many MB each) stay on HPC.
# Only the merged Zea_mays_counts{,_raw}.txt land in the repo, ready to
# git add + push.
#
# Merges every pool_N directory under $baseDir/starsolo/ and stops if any is
# incomplete. To merge a subset, submit with e.g.
#     sbatch --export=ALL,POOLS="1 2 3 4 5 7 8 9" q_04_merge_pools.sh

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq

# STARSOLO_ROOT env var points 04_merge_pools.R at the STARsolo outputs on
# HPC rather than the local (unpopulated) data/starsolo/pool_N tree.
STARSOLO_ROOT="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/starsolo" \
  Rscript scripts/04_merge_pools.R
