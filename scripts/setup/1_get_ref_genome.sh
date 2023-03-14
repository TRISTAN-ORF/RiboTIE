#!/bin/bash
cd ../../data/genome

version="107"

# download primary assembly
mkdir pa
cd pa
wget ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip *.gz

# create sizes.genome file
cd ..
faidx pa/*.fa -i chromsizes > chrom.sizes

# create fa transcriptome
wget http://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz
gunzip *.gz

#gffread -F -w Homo_sapiens.GRCh38.${version}.transcript.fa -g pa/Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.${version}.gtf

# get excl RNA (included in download)
# NCBI search "Homo sapiens"[porgn:__txid9606] AND (biomol_snrna[PROP] OR biomol_snorna[PROP] OR biomol_trna[PROP] OR biomol_rrna[PROP])
