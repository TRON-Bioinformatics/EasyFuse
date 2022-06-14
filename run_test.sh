#!/bin/bash

rm -rf test_out_ini

python processing.py \
    -i test_case/SRR1659960_05pc_R1.fastq.gz \
    test_case/SRR1659960_05pc_R2.fastq.gz \
    -o test_out_ini/ \
    -c config.ini \
    -p 30448
