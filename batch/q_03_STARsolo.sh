#!/bin/bash
#SBATCH --job-name=STARsolo
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=170G
#SBATCH --time=24:00:00
#SBATCH --output=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/STARsolo_pool%a.%A.out
#SBATCH --error=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/STARsolo_pool%a.%A.err

# run from batch/ whatever directory sbatch was called from
cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch

# One STARsolo run per pool. Submit as a job array, one task per pool (pool = plate);
# list exactly the pools to run, from any directory:
#     sbatch --array=1-4 /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/q_03_STARsolo.sh
#     sbatch --array=5,7,8,9 /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/q_03_STARsolo.sh
# SLURM_ARRAY_TASK_ID picks which pool this task processes.
#
# Pool 1 is ~40 GB R1 and takes the longest (~4-6 h). Pools 3-4 finish faster.
# 24 h walltime: the 2026 plates hold ~0.9-1.0 B read pairs each, several
# times an old pool, so the old 8 h limit is not safe for them.

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

: "${SLURM_ARRAY_TASK_ID:?submit with sbatch --array=<pool numbers>}"
bash /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/scripts/03_STARsolo_per_pool.sh "${SLURM_ARRAY_TASK_ID}"
