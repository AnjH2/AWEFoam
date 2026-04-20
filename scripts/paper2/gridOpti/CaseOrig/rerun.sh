#!/bin/bash

#SBATCH --job-name=$(basename "$PWD")
#SBATCH --partition="q64,q36"
#SBATCH --ntasks=2
#SBATCH --mincpus=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time="72:00:00"
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt

# Step 1: Reconstruct the case
reconstructPar

# Step 2: Copy contents of 0 folder to the latest time directory without overwriting
latestTime=$(foamListTimes -case . -latestTime)
cp -rn 0/* $latestTime/

# Step 3: Modify system/controlDict using sed
# Replace startTime with latestTime
timestamp=$(foamListTimes -case . -latestTime)
sed -i "s/^startTime.*/startTime       $timestamp;/" system/controlDict

# Replace maxCo value
sed -i "s/^maxCo[ \t]*[0-9.]*;/maxCo       0.5;/" system/controlDict

# Step 4: Run the final steps
decomposePar -force
srun --resv-ports hisDriftFluxFoam -parallel
reconstructPar

echo "OpenFOAM case rerun complete."
