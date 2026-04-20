#!/bin/bash -l

#SBATCH --job-name="base_rr"
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks=256            # Total number of mpi tasks
#SBATCH --mem=220G
#SBATCH --time="48:00:00"
#SBATCH --partition=standard    # partition name
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --account=project_465002593  # Project for billing
#SBATCH --output=outRe.%j
#SBATCH --error=errRe.%j

ioRanks="($(seq -s ' ' 0 128 "$SLURM_NTASKS"))"
#decomposePar -fileHandler collated -ioRanks "(0 128 $nPoc)" -force


#decomposePar -force -fields
#srun -n $SLURM_NTASKS setFields -parallel
#decomposePar -force -fileHandler collated -ioRanks "$ioRanks"
#srun -n $SLURM_NTASKS topoSet -dict system/topoSetDict.zone -parallel -fileHandler collated -ioRanks "$ioRanks" 
#srun -n $SLURM_NTASKS setFields -parallel -fileHandler collated -ioRanks "$ioRanks" 
echo "simulation started" | tee -a log.ls
srun -n $SLURM_NTASKS alkaWEFoam -parallel -fileHandler collated -ioRanks "$ioRanks" > log.alkaWEFoam
echo "simulation done" | tee -a log.ls

#rm -r pro*


