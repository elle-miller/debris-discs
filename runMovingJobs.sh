#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 10
alpha=4
startingDir=202
n=31
radialRes=200
position=90
for alpha in 3 4; do
  for amplitude in 3 10; do
    for velocity in 0.1 0.3 1 3; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done
