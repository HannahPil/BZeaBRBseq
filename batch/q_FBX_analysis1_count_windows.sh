#!/bin/bash
#SBATCH --job-name=FBX_a1_count
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/FBX_a1_count.%j.out
#SBATCH --error=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/FBX_a1_count.%j.err

# run from batch/ whatever directory sbatch was called from
cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

Rscript /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/scripts/FBX_analysis1_count_windows.R
