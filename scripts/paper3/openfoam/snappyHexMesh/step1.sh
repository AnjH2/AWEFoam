#!/bin/bash -l

#SBATCH --job-name="BASE_step1"
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks=256            # Total number of mpi tasks
#SBATCH --mem=224G
#SBATCH --time="48:00:00"
#SBATCH --partition=standard    # partition name
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --account=project_465002593 # Project for billing
#SBATCH --output=out1.%j
#SBATCH --error=err1.%j


#module load mpi/openmpi/4.0.5-gcc-10.2.0
#source /home/s.hess/OpenFOAM/OpenFOAM-v2106/etc/bashrc

# sbatch createMesh.sh
rm -r 0
cp -rf 0.Orig 0
rm -r result_* log*
#python3 parameter.py

#./createSTL.sh
#python3 genMeshPara.py 1
echo "createMesh" | tee -a log.ls
ls >> log.ls

#makeing MESH
blockMesh
ioRanks="($(seq -s ' ' 0 128 "$SLURM_NTASKS"))"
decomposePar -fileHandler collated -ioRanks "$ioRanks" -force
ls >> log.ls
surfaceFeatureExtract
echo "snappy start" | tee -a log.ls
srun -n $SLURM_NTASKS snappyHexMesh -dict system/snappyHexMeshDict_channel -overwrite -parallel -fileHandler collated -ioRanks "$ioRanks" | tee log.snappyHexMeshDict_channel
echo "snappy done" | tee -a log.ls

#DONE
