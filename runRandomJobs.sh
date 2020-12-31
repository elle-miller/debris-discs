#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Dec 17
startingDir=45
n=30
velocity=0
planOn=1
#alpha=3
#amplitude=2
#width=2
#sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn $invert
#((startingDir+=1))

# First, implement stammlers params and vary position
# A = 3,4 and p = 10,20,30,60,90
invert=0
for alpha in 3 4; do
  for amplitude in 30; do
    for position in 10 20 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn $invert
      ((startingDir+=1))
    done
  done
done

# Then invert the gas bump
# invert gas bump
# p = 50, 100
# A = 3, 10, 30
invert=1
for alpha in 3 4; do
  for amplitude in 3 30; do
    for position in 50 100; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn $invert
      ((startingDir+=1))
    done
  done
done
