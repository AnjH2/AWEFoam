#!/bin/bash -l

#SBATCH --job-name="BASE_step3"
#SBATCH --nodes=2              # Total number of nodes 
#SBATCH --ntasks=256           # Total number of mpi tasks
#SBATCH --mem=224G
#SBATCH --time="48:00:00"
#SBATCH --partition=standard    # partition name
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --account=project_465002593  # Project for billing
#SBATCH --output=out3.%j
#SBATCH --error=err3.%j


#module load mpi/openmpi/4.0.5-gcc-10.2.0
#source /home/s.hess/OpenFOAM/OpenFOAM-v2106/etc/bashrc
ioRanks="($(seq -s ' ' 0 128 "$SLURM_NTASKS"))"

srun -n 1 reconstructParMesh -constant 

checkMesh | tee log.checkMesh1
decomposePar -force -fileHandler collated -ioRanks "$ioRanks"
ls >> log.ls

echo "preprocessing start" | tee -a log.ls

srun -n $SLURM_NTASKS transformPoints -scale '(0.001 0.001 0.001)' -parallel -fileHandler collated -ioRanks "$ioRanks"

srun -n $SLURM_NTASKS topoSet -dict system/topoSetDict.DiaRef -parallel -fileHandler collated -ioRanks "$ioRanks" > log.topoSetDict.DiaRef
srun -n $SLURM_NTASKS refineMesh -overwrite -parallel -fileHandler collated -ioRanks "$ioRanks" &> log.refineMesh

srun -n $SLURM_NTASKS topoSet -dict system/topoSetDict.DiaRefInterFace -parallel -fileHandler collated -ioRanks "$ioRanks" > log.topoSetDict.DiaRefInterFace
srun -n $SLURM_NTASKS refineMesh -dict system/refineMeshDict.DiaRefInterFace -overwrite -parallel -fileHandler collated -ioRanks "$ioRanks" &> log.refineMeshDiaRefInterFace

srun -n $SLURM_NTASKS topoSet -dict system/topoSetDict.zone -parallel -fileHandler collated -ioRanks "$ioRanks" > log.topoSetDict.zone
srun -n $SLURM_NTASKS setFields -parallel -fileHandler collated -ioRanks "$ioRanks" > log.setFields

echo "preprocessing done" | tee -a log.ls
srun -n $SLURM_NTASKS checkMesh -parallel -fileHandler collated -ioRanks "$ioRanks" | tee log.checkMesh2

touch open.foam

#decomposePar -force -fields
echo "simulation started" | tee -a log.ls
srun -n $SLURM_NTASKS alkaWEFoam -parallel -fileHandler collated -ioRanks "$ioRanks" > log.alkaWEFoam
echo "simulation done" | tee -a log.ls
