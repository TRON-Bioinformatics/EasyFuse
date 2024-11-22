#!/bin/bash

cd bin/fusionannotation
python -m venv .venv
source .venv/bin/activate
pip install gffutils
pip install biopython
pip install xxhash

python -m tests.test_breakpoint
python -m tests.test_feature_validation
python -m tests.test_fusion_transcript
python -m tests.test_gff_db_controller
python -m tests.test_io_methods
python -m tests.test_result_handler
python -m tests.test_sequence_handler

python -m src.fusionannotation.fusionannotation \
  --detected_fusions ref_data/Detected_Fusions_miniannotation.csv \
  --annotation_db ref_data/Homo_sapiens.GRCh38.110_miniannotation.gff3.db \
  --tsl_info ref_data/Homo_sapiens.GRCh38.110_miniannotation.gtf.tsl \
  --genome_fasta ref_data/Homo_sapiens.GRCh38.110_minigenome.fa \
  --tsl_filter_level 1,2,NA \
  --out_csv annotation_output.csv
