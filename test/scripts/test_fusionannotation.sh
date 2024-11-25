#!/bin/bash

cd bin
python -m venv fusionannotation/.venv
source fusionannotation/.venv/bin/activate
pip install gffutils
pip install biopython
pip install xxhash

python -m fusionannotation.tests.test_breakpoint
python -m fusionannotation.tests.test_feature_validation
python -m fusionannotation.tests.test_fusion_transcript
python -m fusionannotation.tests.test_gff_db_controller
python -m fusionannotation.tests.test_io_methods
python -m fusionannotation.tests.test_result_handler
python -m fusionannotation.tests.test_sequence_handler

python -m fusionannotator \
  --detected_fusions fusionannotation/ref_data/Detected_Fusions_miniannotation.csv \
  --annotation_db fusionannotation/ref_data/Homo_sapiens.GRCh38.110_miniannotation.gff3.db \
  --tsl_info fusionannotation/ref_data/Homo_sapiens.GRCh38.110_miniannotation.gtf.tsl \
  --genome_fasta fusionannotation/ref_data/Homo_sapiens.GRCh38.110_minigenome.fa \
  --tsl_filter_level 1,2,NA \
  --out_csv fusionannotation/output/annotation_output.csv
