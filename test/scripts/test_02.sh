#!/bin/bash


source test/bin/assert.sh
echo -e "sample_1\t"`pwd`"/test/data/SRR1659960_05pc_R1.fastq.gz\t"`pwd`"/test/data/SRR1659960_05pc_R2.fastq.gz\n" > test/data/test_input.txt
echo -e "sample_2\t"`pwd`"/test/data/SRR1659960_05pc_R1.fastq.gz\t"`pwd`"/test/data/SRR1659960_05pc_R2.fastq.gz" >> test/data/test_input.txt

nextflow main.nf \
  -profile test,conda \
  --output test/output/test2 \
  --input_files test/data/test_input.txt \
  --reference $1 \
  -resume


test -s test/output/test2/sample_1/fusions.csv || { echo "Missing file!"; exit 1; }
test -s test/output/test2/sample_1/fusions.pass.csv || { echo "Missing file!"; exit 1; }
test -s test/output/test2/sample_2/fusions.csv || { echo "Missing file!"; exit 1; }
test -s test/output/test2/sample_2/fusions.pass.csv || { echo "Missing file!"; exit 1; }
