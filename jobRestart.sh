#!/bin/bash -l
#SBATCH --job-name=restart

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

# Shell script to restart lost jobs
# $1 - Directory
# $3 - Amplitude
# $4 - Velocity
# $5 - Position
export OMP_NUM_THREADS=4

srun python ./restart.py -z "$1" -b "$2" -v "$3" -p "$4"
