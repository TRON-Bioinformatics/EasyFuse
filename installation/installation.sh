#!/bin/bash


wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz
wget https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz

tar xvfz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
tar xvfz Homo_sapiens.GRCh38.100.gff3.gz
tar xvfz Homo_sapiens.GRCh38.100.gtf.gz

python easy-fuse gff3-to-db --gff_input Homo_sapiens.GRCh38.100.gff3 --db_output Homo_sapiens.GRCh38.100.gff3.db

python easy-fuse gtf-to-tsl --gtf Homo_sapiens.GRCh38.100.gtf --tsl Homo_sapiens.GRCh38.100.gtf.tsl

