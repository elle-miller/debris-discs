#!/bin/bash -l
#SBATCH --job-name=refined

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

export OMP_NUM_THREADS=4

srun python ./start.py
