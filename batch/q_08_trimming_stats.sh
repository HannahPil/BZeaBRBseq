#!/bin/bash
#BSUB -J 08_trimstats
#BSUB -q sara
#BSUB -n 1
#BSUB -o logs/08_trimstats.%J.out
#BSUB -e logs/08_trimstats.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

../scripts/08_trimming_stats.sh
