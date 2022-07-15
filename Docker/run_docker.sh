#!/bin/bash


# Specify path to reference folder,
REFERENCE_FOLDER=$1

# Specify path to data folder,
DATA_FOLDER=$2

# Specify image tag, e.g. `docker.io/tronbioinformatics/easyfuse:latest`
IMAGE_TAG=$3

# Specify output folder name, e.g. `results`
OUTPUT_FOLDER_NAME=$4

docker run \
--name easyfuse_container \
-v ${REFERENCE_FOLDER}:/ref \
-v ${DATA_FOLDER}:/data \
--rm \
-it ${IMAGE_TAG} \
python /code/easyfuse/processing.py -i /data/input_fastqs -o /data/${OUTPUT_FOLDER}/
