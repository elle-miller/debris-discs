#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: June 7

for z in 107 119 123 41 43 44 45 46 47 48 49; do
  sbatch jobRestart.sh $z
done