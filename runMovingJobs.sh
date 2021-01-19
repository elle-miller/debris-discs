#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Jan 19
alpha=4
startingDir=130
n=31
radialRes=100
position=90
for alpha in 3; do
  for amplitude in 3 30; do
    for velocity in 0.1 0.3 1 3; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done
for alpha in 4; do
  for amplitude in 3 30; do
    for velocity in 1 3 10 30; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done
