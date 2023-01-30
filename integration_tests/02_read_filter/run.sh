#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/02_read_filter
input=$test_folder/input.SRR1659960_05pc_Aligned.out.downsampled.bam
output=$test_folder/observed.SRR1659960_05pc_Aligned.out.downsampled.bam

easy-fuse read-filter \
--input $input \
--output $output

test -s $output || { echo "Missing BAM output file!"; exit 1; }
assert_eq `samtools view $output | wc -l` 4600 "Wrong number of reads"
