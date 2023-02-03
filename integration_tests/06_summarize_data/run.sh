#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/06_summarize_data
input_fusions=$test_folder/input.fusions.csv
input_fusions_contexts=$test_folder/input.fusions_context_sequences.csv
input_requant_cpm=$test_folder/input.requant_cpm.tdt
input_requant_counts=$test_folder/input.requant_counts.counts
input_reads_stats=$test_folder/input.reads_stats.txt

mkdir $test_folder/observed

easy-fuse summarize-data \
  --input-fusions $input_fusions \
  --input-fusion-context-seqs $input_fusions_contexts \
  --input-requant-cpm $input_requant_cpm \
  --input-requant-counts $input_requant_counts \
  --input-reads-stats $input_reads_stats \
  --model_predictions \
  --output-folder $test_folder/observed \
  --requant-mode best \
  --model-pred-threshold 0.5

test -s $test_folder/observed/fusions.pass.csv || { echo "Missing predicted fusions output file!"; exit 1; }
cat $test_folder/observed/fusions.pass.csv
assert_eq `wc -l $test_folder/observed/fusions.pass.csv | cut -d " " -f 1` 170 "Wrong number of fusions"
#
