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

#Running case
hisDriftFluxFoam

reconstructPar
echo "Setup and initial run steps completed."

