#!/bin/bash


#source bin/assert.sh
output=output/test1

echo -e "sample_name\t"`pwd`"/test/data/SRR1659960_05pc_R1.fastq.gz\t"`pwd`"/test/data/SRR1659960_05pc_R2.fastq.gz" > test/data/test_input.txt
nextflow main.nf \
  -profile test,conda \
  --reference /projects/data/human/ensembl/GRCh38.86/STAR_idx/ \
  --output $output \
  --input_files test/data/test_input.txt

#test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }
