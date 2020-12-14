#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Dec 14
startingDir=30
n=100
velocity=0
for alpha in 3 4; do
  for amplitude in 10 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn
      ((startingDir+=1))
    done
  done
done
