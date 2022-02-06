#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 25 vfrag for 1e-3

startingDir=240
n=31
radialRes=200
alpha=3
velocity=0
deltaRZ=0 # multiplying factor not 0, but 10^0 = 1
deltaT=0
for amplitude in 10 30; do
  for position in 30 60 90; do
    for vfrag in 100 300; do
      sbatch jobFrag.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $vfrag $deltaRZ $deltaT
      ((startingDir+=1))
    done
  done
done

# Now run 189 case with reduced delta for stationary
amplitude=10
velocity=0
position=90
vfrag=100
deltaT=2 # mulitplying factor 10^-2 to give 10^-5
for deltaR in 0 2; do
  sbatch jobFrag.sh "$startingDir" $alpha $amplitude $velocity $position $n $radialRes $vfrag $deltaRZ $deltaT
  ((startingDir+=1))
done
