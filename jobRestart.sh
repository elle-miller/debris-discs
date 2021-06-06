#!/bin/bash -l
#SBATCH --job-name=restart

module load anaconda3
module load openmpi-4.0.1

# Shell script to restart lost jobs
export OMP_NUM_THREADS=4

srun python ./restart.py -z "$1"
