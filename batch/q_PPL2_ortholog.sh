#!/bin/bash
#BSUB -J PPL2_orth
#BSUB -q sara
#BSUB -n 4
#BSUB -o logs/PPL2_orth.%J.out
#BSUB -e logs/PPL2_orth.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00

# blastp / makeblastdb come from the cluster module, not the pipeline conda env
# (the maize env has no BLAST). seqtk is absent too; the script falls back to
# awk for sequence extraction, so nothing else is needed here.
module load blast/2.17.0

bash ../scripts/PPL2_ortholog_check_hpc.sh
