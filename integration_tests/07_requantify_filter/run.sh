#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/07_requantify_filter
input=$test_folder/input.SRR1659960_05pc_Aligned.out.downsampled.bam
input2=$test_folder/input.context_seqs.csv.debug
output=$test_folder/observed.bam

easy-fuse requantify-filter \
  --input $input \
  --input2 $input2 \
  --input-reads-stats $test_folder/input.reads_stats.txt \
  --output $output

test -s $output || { echo "Missing filtered BAM output file!"; exit 1; }
assert_eq `samtools view $output | wc -l` 3404 "Wrong number of reads"
