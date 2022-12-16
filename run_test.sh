#!/bin/bash

rm -rf test_out_slurm

python processing.py \
    -i test_case/SRR1659960_05pc_R1.fastq.gz \
    test_case/SRR1659960_05pc_R2.fastq.gz \
    -o test_out_slurm/ \
    -c config.ini \
    -s 30448 \
    -q slurm \
    -p Compute


rm -rf test_out_shell

python processing.py \
    -i test_case/SRR1659960_05pc_R1.fastq.gz \
    test_case/SRR1659960_05pc_R2.fastq.gz \
    -o test_out_shell/ \
    -c config.ini \
    -s 30448 \
    -q none
