#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Jan 21 ~ high res with fixed PF (less sharp with tanh, cancelled 5 last bins)

startingDir=150
n=31
velocity=0
radialRes=300
massRes=200
for alpha in 3 4; do
  for amplitude in 3 10 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $massRes
      ((startingDir+=1))
    done
  done
done

