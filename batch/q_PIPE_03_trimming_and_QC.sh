#!/bin/bash
#BSUB -J 03_trimming
#BSUB -q sara
#BSUB -n 1
#BSUB -o logs/03_trimming.%J.out
#BSUB -e logs/03_trimming.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00


../scripts/PIPE_03_trimming_and_QC.sh
