#!/bin/bash

python -m venv bin/fusionannotation/.venv
source bin/fusionannotation/.venv/bin/activate
pip install gffutils
pip install biopython
pip install xxhash

python -m unittest bin/fusionannotation/tests/test_breakpoint.py
python -m unittest bin/fusionannotation/tests/test_feature_validation.py
python -m unittest bin/fusionannotation/tests/test_fusion_transcript.py
python -m unittest bin/fusionannotation/tests/test_gff_db_controller.py
python -m unittest bin/fusionannotation/tests/test_io_methods.py
python -m unittest bin/fusionannotation/tests/test_result_handler.py
python -m unittest bin/fusionannotation/tests/test_sequence_handler.py
