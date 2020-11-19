#!/bin/bash -l
#SBATCH --job-name=loop
#SBATCH --mail-user elle.miller101@gmail.com
#SBATCH --mail-type=ALL

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

# Shell script to run individual jobs. Current inputs:
# $1 - Directory
# $2 - Alpha
# $3 - Amplitude
# $4 - Velocity
# $5 - Position

srun python ./main.py -z "$1" -a 1e-"$2" -b "$3" -v "$4" -p "$5" -n "$6" -4 "$7"
