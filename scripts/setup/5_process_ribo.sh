#!/bin/bash

cd ../../data/

h5='GRCh38v107.h5'
ribo_meta='ribo/metadata.txt'

python ../scripts/setup/s_process_ribo_reads.py $h5 $ribo_meta