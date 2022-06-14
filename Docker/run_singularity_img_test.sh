#!/bin/bash


singularity exec docker://docker.io/tronbioinformatics/easyfuse:1.3.5 \
--bind /scratch/info/data/easyfuse/easyfuse_ref:/ref \
--bind /scratch/info/data/easyfuse/easyfuse_data:/data \
python /code/easyfuse/processing.py -i /data/input_fastqs -o /data/results/

