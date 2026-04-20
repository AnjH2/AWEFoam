#!/bin/bash -l
#SBATCH --job-name=BASE_step2  # Job name
#SBATCH --output=out2.%j # Name of stdout output file
#SBATCH --error=err2.%j  # Name of stderr error file
#SBATCH --partition=small       # Partition name
#SBATCH --ntasks=1              # One task (process)
#SBATCH --nodes=1
#SBATCH --mem=224G
#SBATCH --time=48:00:00         # Run time (hh:mm:ss)
#SBATCH --account=project_465002593  # Project for billing


ioRanks="($(seq -s ' ' 0 128 "128"))"
srun -n $SLURM_NTASKS reconstructParMesh -constant 

checkMesh | tee log.checkMesh1
decomposePar -force -fileHandler collated -ioRanks "$ioRanks"
ls >> log.ls

