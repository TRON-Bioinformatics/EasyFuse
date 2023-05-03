#!/bin/bash


source test/bin/assert.sh
echo -e "sample_name\t"`pwd`"/test/data/SRR1659960_05pc_R1.fastq.gz\t"`pwd`"/test/data/SRR1659960_05pc_R2.fastq.gz" > test/data/test_input.txt

nextflow main.nf \
  -profile test,conda \
  --output test/output/test1 \
  --input_files test/data/test_input.txt \
  --reference /scratch/info/data/easyfuse/easyfuse_ref/ \
  -resume


test -s test/output/test1/sample_name/fusions_1.csv || { echo "Missing file!"; exit 1; }
test -s test/output/test1/sample_name/fusions.pass.csv || { echo "Missing file!"; exit 1; }
test -s test/output/test1/sample_name/fusions.pass.all.csv || { echo "Missing file!"; exit 1; }
