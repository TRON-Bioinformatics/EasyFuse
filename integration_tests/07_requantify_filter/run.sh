#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/07_requantify_filter
input=integration_tests/02_read_filter/input.SRR1659960_05pc_Aligned.out.downsampled.bam
input2=$test_folder/input.reads_stats.txt
output=$test_folder/observed.bam

easy-fuse requantify-filter \
  --input $input \
  --input2 $input2 \
  --output $output

test -s $output || { echo "Missing filtered BAM output file!"; exit 1; }
assert_eq `samtools view $output | wc -l` 2 "Wrong number of reads"
