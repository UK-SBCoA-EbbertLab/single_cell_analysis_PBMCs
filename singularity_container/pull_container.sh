#!/bin/bash

# analysis container
#singularity pull --arch amd64 

# Pull the latest version of the single-cell analysis container
singularity pull library://phdo222/single_cell_data/single_cell_analysis:2025-03-04
