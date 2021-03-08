#!/bin/bash -l
#SBATCH --job-name=invert

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

# Shell script to run individual jobs. Current inputs:
# $1 - Directory
# $2 - Alpha
# $3 - Amplitude
# $4 - Velocity
# $5 - Position
# $6 - Number of snapshots
# $7 - Radial resolution
# $8 - Invert
export OMP_NUM_THREADS=4


srun python ./main.py -z "$1" -a 1e-"$2" -b "$3" -v "$4" -p "$5" -n "$6" -r "$7" -i "$8"