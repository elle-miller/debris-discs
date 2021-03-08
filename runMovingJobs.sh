#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Sharp formation, tanh commented out
startingDir=270
n=31
radialRes=200
# position=90
amplitude=10
velocity=1
alpha=3
for position in 60 80 85; do
  sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
  ((startingDir+=1))
done

#for alpha in 3 4; do
#  for amplitude in 3 10; do
#    for velocity in 0.1 0.3 1 3; do
#      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#      ((startingDir+=1))
#    done
#  done
#done

# Main set
#for alpha in 3 4; do
#  for amplitude in 3 10; do
#    for velocity in 0.1 0.3 1 3; do
#      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#      ((startingDir+=1))
#    done
#  done
#done
