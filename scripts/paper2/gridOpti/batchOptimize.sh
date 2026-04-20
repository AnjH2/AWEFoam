#!/bin/bash

#SBATCH --job-name="Search Manager G"
#SBATCH --partition="q64,q36"
#SBATCH --ntasks=1
#SBATCH --mincpus=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --time="168:00:00"
#SBATCH --mail-user="anj@mpe.au.dk"
#SBATCH --mail-type=all
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt


python3 OManagerStuktured.py 
