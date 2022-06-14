#!/bin/bash

docker run \
--name test_easyfuse_container \
-v /scratch/info/data/easyfuse/easyfuse_ref:/ref \
-v /scratch/info/data/easyfuse/easyfuse_data:/data \
--rm \
-it docker.io/tronbioinformatics/easyfuse:1.3.5 \
python /code/easyfuse/processing.py -i /data/input_fastqs -o /data/results/
