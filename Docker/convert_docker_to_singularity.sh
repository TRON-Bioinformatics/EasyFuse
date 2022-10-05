#!/bin/bash

# Specify image ID from `docker images`
IMAGE_ID=$1

# Specify output folder to store files
OUTPUT_DIR=$2

# Save Docker image to .tar archive
docker save <image_id> -o ${OUTPUT_DIR}/easyfuse_latest.tar

# Build Singularity image from .tar archive and save it to .sif file
singularity build ${OUTPUT_DIR}/easyfuse_latest.sif docker-archive://${OUTPUT_DIR}/easyfuse_latest.tar
