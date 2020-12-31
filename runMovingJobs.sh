#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Nov 8 - running on Nov19
alpha=4
startingDir=43

# Loop through two alphas, two amplitudes and three positions
n=50
position=90
for alpha in 3 4; do
  for amplitude in 10 30; do
    for velocity in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn
      ((startingDir+=1))
    done
  done
done
