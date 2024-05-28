#!/bin/bash

# This script sets up and runs an OpenFOAM case

# Echo commands to terminal
set -x

# Copy the initial conditions
cp -r 0.orig 0

# Generate the mesh
blockMesh

# Set the initial field values
setFields

# Decompose the domain for parallel run
decomposePar -force

#Running case
mpirun -np 4 hisDriftFluxFoam -parallel

echo "Setup and initial run steps completed."

