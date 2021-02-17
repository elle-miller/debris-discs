#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 17 rerun 189/208 with lots more Nt

startingDir=235
n=300
radialRes=200
alpha=3
amplitude=10
position=90
velocity=0
for v in 0 1; do
  sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
  ((startingDir+=1))
done

#for amplitude in 10 30; do
#  for position in 30 60 90; do
#    for vfrag in 3 10; do
#      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $vfrag
#      ((startingDir+=1))
#    done
#  done
#done
#for alpha in 3 4; do
#  for amplitude in 3 10 30; do
#    for position in 30 60 90; do
#      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#      ((startingDir+=1))
#    done
#  done
#done



