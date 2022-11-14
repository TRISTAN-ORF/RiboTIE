#!/bin/bash

cd ../../data/genome
version="107"
mkdir STAR

# setting up reference genome for STAR mapping
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir STAR/ --genomeFastaFiles pa/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.${version}.gtf --sjdbOverhang 35

# set up redundant RNA folder
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir STAR/excl_RNA --genomeFastaFiles excl_RNA.fasta --genomeSAindexNbases 6