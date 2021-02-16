#!/bin/bash -l
#SBATCH --job-name=widebelt

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

# Shell script to run individual jobs. Current inputs:
# $1 - Directory
# $2 - Alpha
# $3 - Amplitude
# $4 - Velocity
# $5 - Position
# $6 - TimeBumpMove
# $7 - Start year power
export OMP_NUM_THREADS=4

srun python ./main.py -z "$1" -a 1e-"$2" -b "$3" -v "$4" -p "$5" -t "$6" -s "$7"
