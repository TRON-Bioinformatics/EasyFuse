#!/bin/bash

# Create installation folder
mkdir -p installation/

# Download and extract references
if [ ! -f installation/Homo_sapiens.GRCh38.dna.primary_assembly.fa ]; then
    wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P installation/
    gunzip installation/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi

if [ ! -f installation/Homo_sapiens.GRCh38.100.gff3 ]; then
    wget https://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz -P installation/
    gunzip installation/Homo_sapiens.GRCh38.100.gff3.gz
fi

if [ ! -f installation/Homo_sapiens.GRCh38.100.gtf ]; then
    wget https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz -P installation/
    gunzip installation/Homo_sapiens.GRCh38.100.gtf.gz
fi

# for starfusion
if [ ! -f installation/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play ]; then
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz -P installation/
    tar xvfz installation/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
fi

# for fusioncatcher
if [ ! -d installation/fusioncatcher_index/ ]; then
    wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa -P installation/
    wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab -P installation/
    wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac -P installation/
    wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad -P installation/
    wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.md5 -P installation/
    md5sum -c installation/human_v102.md5
    if [ "$?" -ne "0" ]; then
	    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files from above have errors! MD5 checksums do not match! Please, download them again or re-run this script again!   \033[0m\n"
	    exit 1
    fi
    cat installation/human_v102.tar.gz.* > installation/human_v102.tar.gz
    rm -f installation/human_v102.tar.gz.*
    mkdir -p installation/fusioncatcher_index/
    if ! tar -xzf installation/human_v102.tar.gz -C .; then
	    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files are corrupted! Please, download them again or re-run this script again!   \033[0m\n"
	    exit 1
    fi
    rm -f installation/human_v102.tar.gz
    rm -f installation/human_v102.md5
fi

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

