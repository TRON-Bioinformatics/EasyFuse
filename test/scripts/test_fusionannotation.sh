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
