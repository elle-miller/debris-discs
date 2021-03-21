#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Sharp formation, tanh commented out
startingDir=290
n=2
radialRes=200
position=90
amplitude=10
velocity=1
alpha=3
for d2g in 0.003 0.03; do
  for zeta in 0.3 0.1 0.01 0.001 0.0001 0.00001; do
    sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $d2g $zeta
    ((startingDir+=1))
  done
done

