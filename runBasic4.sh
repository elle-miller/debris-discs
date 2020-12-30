#!/bin/bash -l
#SBATCH --job-name=basic4
#SBATCH --mail-user elle.miller101@gmail.com
#SBATCH --mail-type=ALL

# Allocate 4 GB of Ram per node
# --mem-per-cpu can be used to allocate RAM per CPU
#SBATCH --mem=4000

# Run 4 tasks (copies) of this job, use this for parallel jobs
# The default is one task per node, can be overwritten with --cpus-per-task
#SBATCH -n 4

module load anaconda3
module load openmpi-4.0.1
source activate dustpy2

srun python ./mainBasic.py -z 87
