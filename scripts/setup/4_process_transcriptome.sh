#!/bin/bash

cd ../../data/
mkdir genome/riboformer

gtf='genome/Homo_sapiens.GRCh38.107.gtf'
fasta='genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
chromsizes='genome/chrom.sizes'
h5_out='GRCh38_v107'
meta_dir='genome/riboformer'


python ../scripts/setup/s_process_transcriptome.py $gtf $fasta $chromsizes $h5_out $meta_dir
