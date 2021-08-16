#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Aug 16 - FINAL JOBS!!!!
startingDir=500
n=31
radialRes=240
amplitude=10
alpha=3
for position in 120 180; do
  for velocity in 0 0.1 0.3 1; do
    sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
    ((startingDir+=1))
  done
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
