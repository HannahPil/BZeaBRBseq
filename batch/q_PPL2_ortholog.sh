#!/bin/bash
#SBATCH --job-name=PPL2_orth
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/PPL2_orth.%j.out
#SBATCH --error=/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch/logs/PPL2_orth.%j.err

# run from batch/ whatever directory sbatch was called from
cd /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/batch

# blastp / makeblastdb come from the cluster module, not the pipeline conda env
# (the maize env has no BLAST). seqtk is absent too; the script falls back to
# awk for sequence extraction, so nothing else is needed here.
module load blast/2.17.0

bash /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/BZeaBRBseq/scripts/PPL2_ortholog_check_hpc.sh
