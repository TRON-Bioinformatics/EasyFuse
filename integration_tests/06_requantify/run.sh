#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/06_requantify
input=integration_tests/02_read_filter/input.SRR1659960_05pc_Aligned.out.downsampled.bam
input_reads_stats=$test_folder/input.reads_stats.txt
output=$test_folder/observed.requantified.tdt

easy-fuse requantify \
  -i $input \
  -o $output \
  -d 10 \
  --input-reads-stats $input_reads_stats

test -s $output || { echo "Missing requantified output file!"; exit 1; }
assert_eq `wc -l` 100 "Wrong number of reads"
