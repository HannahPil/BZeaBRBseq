#!/bin/bash
#BSUB -J REC_demux
#BSUB -q sara
#BSUB -n 4
#BSUB -o REC_demux.%J.out
#BSUB -e REC_demux.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -W 6:00

../scripts/REC_recover_missing_demux.sh
