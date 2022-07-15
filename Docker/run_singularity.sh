#!/bin/bash


# Specify path to reference folder,
REFERENCE_FOLDER=$1

# Specify path to data folder,
DATA_FOLDER=$2

# Specify image tag or SIF, e.g. `docker://docker.io/tronbioinformatics/easyfuse:latest`
DOCKERHUB_OR_SIF=$3

# Specify output folder name, e.g. `results`
OUTPUT_FOLDER_NAME=$4


singularity exec \
--containall \
--bind ${REFERENCE_FOLDER}:/ref \
--bind ${DATA_FOLDER}:/data \
${DOCKERHUB_OR_SIF} \
python /code/easyfuse/processing.py -i /data/input_fastqs -o /data/${OUTPUT_FOLDER_NAME}/
