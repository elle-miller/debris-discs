#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 1 ~ Nr=200, 0.03 denom + mmax=1e10

startingDir=184
n=31
velocity=0
radialRes=200
for alpha in 3 4; do
  for amplitude in 3 10 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done

