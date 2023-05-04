#!/bin/bash

# Create installation folder
mkdir -p installation/

# Download references
wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P installation/
wget https://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz -P installation/
wget https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz -P installation/

# Extract references
gunzip installation/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip installation/Homo_sapiens.GRCh38.100.gff3.gz
gunzip installation/Homo_sapiens.GRCh38.100.gtf.gz

# Convert GFF3 to DB format
easy-fuse gff3-to-db --gff_input installation/Homo_sapiens.GRCh38.100.gff3 --db_output installation/Homo_sapiens.GRCh38.100.gff3.db

# Convert GTF to TSL format
easy-fuse gtf-to-tsl --gtf installation/Homo_sapiens.GRCh38.100.gtf --tsl installation/Homo_sapiens.GRCh38.100.gtf.tsl

# Generate STAR index
STAR \
    --runMode genomeGenerate \
    --runThreadN 12 \
    --genomeDir installation/star_index \
    --genomeFastaFiles installation/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile installation/Homo_sapiens.GRCh38.100.gtf

# Generate Fusioncatcher reference
# fusioncatcher-build ...

# Generate Starfusion reference
# STAR-Fusion ...
