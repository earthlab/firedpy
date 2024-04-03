#!/bin/bash
# entrypoint.sh

# Activate the conda environment
source activate firedpy

# Execute the command provided to the docker run command
exec "$@"
