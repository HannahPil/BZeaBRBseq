#!/bin/bash
#BSUB -J 02_demux
#BSUB -q sara
#BSUB -n 4
#BSUB -o 02_demux.%J.out
#BSUB -e 02_demux.%J.err
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16000]"
#BSUB -W 24:00

../scripts/02_demultiplex.sh
