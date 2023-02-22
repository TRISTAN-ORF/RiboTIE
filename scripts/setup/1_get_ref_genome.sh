#!/bin/bash
cd ../../data/genome

version="107"

# download primary assembly
wget ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip *.gz

# create sizes.genome file
faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa -i chromsizes > chrom.sizes

# create fa transcriptome
wget http://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz
gunzip ftp://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz

# get excl RNA (included in download)
# NCBI search "Homo sapiens"[porgn:__txid9606] AND (biomol_snrna[PROP] OR biomol_snorna[PROP] OR biomol_trna[PROP] OR biomol_rrna[PROP])