#!/bin/bash
#BSUB -J 03_trimming
#BSUB -q sara
#BSUB -n 1
#BSUB -o 03_trimming.%J.out
#BSUB -e 03_trimming.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00


../scripts/03_trimming_and_QC.sh
