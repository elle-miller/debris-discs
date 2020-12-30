#!/bin/bash -l
#SBATCH --job-name=basic1
#SBATCH --mail-user elle.miller101@gmail.com
#SBATCH --mail-type=ALL

# Nothing new - control test

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

srun python ./mainBasic.py -z 84
