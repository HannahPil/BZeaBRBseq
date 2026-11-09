#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/trim_pool%a.%A.out
#SBATCH --error=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/trim_pool%a.%A.err

# run from batch/ whatever directory sbatch was called from
cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch

# Submit as a job array, one task per pool (pool = plate);
# list exactly the pools to run, from any directory:
#     sbatch --array=1-4 /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/q_02b_trim.sh
#     sbatch --array=5,7,8,9 /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/q_02b_trim.sh
# Pool 1 takes longest of plates 1-4 (~40 GB R1); the 2026 plates are larger,
# hence the 24 h walltime.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

: "${SLURM_ARRAY_TASK_ID:?submit with sbatch --array=<pool numbers>}"
bash /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/scripts/02b_trim_pools.sh "${SLURM_ARRAY_TASK_ID}"
