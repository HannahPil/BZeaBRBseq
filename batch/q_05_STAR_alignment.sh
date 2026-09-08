#!/bin/bash
#BSUB -J 05_STAR
#BSUB -q sara
#BSUB -n 12
#BSUB -o 05_STAR.%J.out
#BSUB -e 05_STAR.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 18:00

../scripts/05_STAR_alignment_GEMMA.sh Zea_mays
