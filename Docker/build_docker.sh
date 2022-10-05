#!/bin/bash

# Specify image tag, e.g. `docker.io/tronbioinformatics/easyfuse:latest`
IMAGE_TAG=$1

# Specify path containing Dockerfile and context
PATH=$2

docker build \
--rm \
-t ${IMAGE_TAG} \
${PATH}
