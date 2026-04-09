#!/bin/bash
#BSUB -J 05b_STAR_cp
#BSUB -q sara
#BSUB -n 12
#BSUB -o 05b_STAR_cp.%J.out
#BSUB -e 05b_STAR_cp.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 6:00

../scripts/05b_STAR_alignment_chloroplast.sh
