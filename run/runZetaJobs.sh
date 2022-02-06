#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Sharp formation, tanh commented out
startingDir=50
n=31
startyear=6
radialRes=240

# Prime stuff
position=90
amplitude=10
velocity=1
alpha=3


steep=0.03
for zeta in 0.3 0.1 0.01 0.001 0.0001 0.00001; do
  sbatch jobZeta.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $steep $zeta $startyear
  ((startingDir+=1))
done

zeta=0.1
for steep in 0.01 0.1; do
  sbatch jobZeta.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $steep $zeta $startyear
  ((startingDir+=1))
done

