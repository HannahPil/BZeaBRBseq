#!/bin/bash
#BSUB -J 07_summary
#BSUB -q sara
#BSUB -n 1
#BSUB -o 07_summary.%J.out
#BSUB -e 07_summary.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00

../scripts/07_generate_summary_statistics.sh
