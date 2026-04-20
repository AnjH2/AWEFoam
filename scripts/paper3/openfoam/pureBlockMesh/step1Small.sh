#!/bin/bash -l

#SBATCH --job-name="100_3_12_run"
#SBATCH --ntasks=46            # Total number of mpi tasks
#SBATCH --mem=79G
#SBATCH --time="48:00:00"
#SBATCH --partition=small    # partition name
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --account=project_465002593 # Project for billing
#SBATCH --output=out1.%j
#SBATCH --error=err1.%j


ioRanks="($(seq -s ' ' 0 128 "$SLURM_NTASKS"))"
python3 makeGrading.py
rm -r 0
cp -rf 0.Orig 0
blockMesh
topoSet -dict system/topoSetDict.outlets
refineMesh -overwrite
#topoSet -dict system/topoSetDict2
#refineMesh -overwrite 
#checkMesh -latestTime -writeAllFields
renumberMesh -overwrite
topoSet -dict system/topoSetDict.zone
setFields
transformPoints -rotate-z -90
checkMesh -latestTime > log.checkMesh

decomposePar -force -fileHandler collated -ioRanks "$ioRanks"
srun -n $SLURM_NTASKS alkaWEFoam -parallel -fileHandler collated -ioRanks "$ioRanks" > log.alkaWEFoam

