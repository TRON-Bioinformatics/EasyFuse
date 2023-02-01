#!/bin/bash


source integration_tests/assert.sh
test_folder=integration_tests/03_fusion_parser
input_fusioncatcher=$test_folder/input_fusioncatcher.summary_candidate_fusions.txt
input_fusioncatcher2=$test_folder/input_fusioncatcher.final-list_candidate-fusion-genes.txt
input_infusion=$test_folder/input_infusion.fusions.detailed.txt
input_mapsplice=$test_folder/input_mapsplice.fusions_well_annotated.txt
input_soapfuse=$test_folder/input_soapfuse.final.Fusion.specific.for.genes
input_starfusion=$test_folder/input_starfusion.star-fusion.fusion_predictions.tsv

output_folder=$test_folder/observed
mkdir $output_folder

output_file=${output_folder}/Detected_Fusions.csv
output_file_fusioncatcher=${output_folder}/fusioncatcher.csv
output_file_infusion=${output_folder}/infusion.csv
output_file_mapsplice=${output_folder}/mapsplice.csv
output_file_soapfuse=${output_folder}/soapfuse.csv
output_file_starfusion=${output_folder}/starfusion.csv

easy-fuse fusion-parser \
--input-fusioncatcher $input_fusioncatcher \
--input-fusioncatcher2 $input_fusioncatcher2 \
--input-infusion $input_infusion \
--input-mapsplice $input_mapsplice \
--input-soapfuse $input_soapfuse \
--input-starfusion $input_starfusion \
--output $output_folder \
--sample test

test -s $output_file || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file | cut -d " " -f 1` 32 "Wrong number of fusions"

test -s $output_file_fusioncatcher || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file_fusioncatcher | cut -d " " -f 1` 9 "Wrong number of fusions"

test -s $output_file_infusion || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file_infusion | cut -d " " -f 1` 1 "Wrong number of fusions"

test -s $output_file_mapsplice || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file_mapsplice | cut -d " " -f 1` 8 "Wrong number of fusions"

test -s $output_file_soapfuse || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file_soapfuse | cut -d " " -f 1` 9 "Wrong number of fusions"

test -s $output_file_starfusion || { echo "Missing CSV output file!"; exit 1; }
assert_eq `wc -l $output_file_starfusion | cut -d " " -f 1` 9 "Wrong number of fusions"
