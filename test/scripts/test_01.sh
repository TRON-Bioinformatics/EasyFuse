#!/bin/bash

echo -e "sample_name\t"`pwd`"/test/data/SRR1659960_05pc_R1.fastq.gz\t"`pwd`"/test/data/SRR1659960_05pc_R2.fastq.gz" > test/data/test_input.txt

nextflow main.nf \
  -profile test,conda \
  --output test/output/test1 \
  --input_files test/data/test_input.txt \
  -resume
