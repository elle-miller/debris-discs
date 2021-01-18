#!/bin/bash -l
#SBATCH --job-name=testPlanOff
#SBATCH --mail-user elle.miller101@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -n 1

export OMP_NUM_THREADS=4

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

srun python ./main.py -n 31 -a 1e-4 -b 30 -p 30 -r 300 -e 6 -4 0 -z 127
