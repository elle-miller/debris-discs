#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Jan 12

startingDir=112
n=31
velocity=0
radialRes=300
for alpha in 3; do
  for amplitude in 3 10 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done
for alpha in 4; do
  for amplitude in 3 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done

