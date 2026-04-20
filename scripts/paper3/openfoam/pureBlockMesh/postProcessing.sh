#!/bin/bash -l
#SBATCH --job-name="BASE_PP"
#SBATCH --ntasks=46           # Total number of mpi tasks
#SBATCH --mem=224G
#SBATCH --time="48:00:00"
#SBATCH --partition=small    # partition name
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --account=project_465002593  # Project for billing


ioRanks="($(seq -s ' ' 0 128 "$SLURM_NTASKS"))"
#srun -n $SLURM_NTASKS postProcess -func extractLine -latestTime -parallel -fileHandler collated -ioRanks "$ioRanks"
srun -n $SLURM_NTASKS postProcess -func surfaces -latestTime -parallel -fileHandler collated -ioRanks "$ioRanks"
srun -n $SLURM_NTASKS foamToVTK -parallel -no-point-data -latestTime
