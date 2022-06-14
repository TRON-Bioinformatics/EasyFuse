#!/bin/bash

docker run \
--name test_easyfuse_container \
-v /scratch/info/data/easyfuse/easyfuse_ref:/ref \
-it docker.io/tronbioinformatics/easyfuse:1.3.5 \
/bin/bash
