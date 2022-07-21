#!/bin/bash

# Specify path to reference folder,
REFERENCE_FOLDER=$1

# Specify image tag, e.g. `docker.io/tronbioinformatics/easyfuse:latest`
IMAGE_TAG=$2


docker run \
--name easyfuse_container \
-v ${REFERENCE_FOLDER}:/ref \
-it ${IMAGE_TAG} \
/bin/bash
