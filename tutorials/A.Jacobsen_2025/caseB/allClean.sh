#!/bin/bash

# This script cleans up an OpenFOAM case

# Echo commands to terminal
set -x

# Remove the 0 directory
rm -rf 0

# Clean mesh files
rm -rf constant/polyMesh
rm -f constant/polyMesh/*


# Optionally, clean any solver-specific directories or files
rm -rf postProcessing*

rm -rf dynamicCode*

echo "Cleanup completed."

