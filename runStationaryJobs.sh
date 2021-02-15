#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 16 - different vfrag parameters

startingDir=220
n=31
velocity=0
radialRes=200
alpha=4
for amplitude in 10 30; do
  for position in 30 60 90; do
    for vfrag in 3 10; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $vfrag
      ((startingDir+=1))
    done
  done
done
#for alpha in 3 4; do
#  for amplitude in 3 10 30; do
#    for position in 30 60 90; do
#      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#      ((startingDir+=1))
#    done
#  done
#done



