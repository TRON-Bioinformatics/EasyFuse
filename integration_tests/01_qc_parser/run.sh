#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/01_qc_parser
input_r1=$test_folder/input.fastqc_data_R1.txt
input_r2=$test_folder/input.fastqc_data_R2.txt

easy-fuse qc-parser \
-i $input_r1 $input_r2 \
-o $test_folder/observed.qc_table.txt

test -s $test_folder/observed.qc_table.txt || { echo "Missing QC table output file!"; exit 1; }
cmp --silent <(cut -f1 -d"," --complement  $test_folder/expected.qc_table.txt) <(cut -f1 -d"," --complement $test_folder/observed.qc_table.txt)
