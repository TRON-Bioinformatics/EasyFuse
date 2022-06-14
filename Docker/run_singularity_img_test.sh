#!/bin/bash


singularity exec \
-W ./ \
--bind /scratch/info/data/easyfuse/easyfuse_ref:/ref \
--bind /scratch/info/data/easyfuse/easyfuse_data:/data \
docker://docker.io/tronbioinformatics/easyfuse:latest \
python /code/easyfuse/processing.py -i /data/input_fastqs -o /data/results/

