#!/bin/bash

python -m unittest bin/fusionannotation/tests/test_breakpoint.py
python -m unittest bin/fusionannotation/tests/test_feature_validation.py
python -m unittest bin/fusionannotation/tests/test_fusion_transcript.py
python -m unittest bin/fusionannotation/tests/test_gff_db_controller.py
python -m unittest bin/fusionannotation/tests/test_io_methods.py

