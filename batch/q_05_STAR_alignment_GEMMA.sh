#!/bin/bash
#BSUB -J 05_STAR
#BSUB -q sara
#BSUB -n 12
#BSUB -o 05_STAR.%J.out
#BSUB -e 05_STAR.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

../scripts/05_STAR_alignment.sh Zea_mays
