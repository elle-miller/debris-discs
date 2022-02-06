#!/bin/bash -l

# FINAL RUNTHROUGH OF EVERYTHING!!!

# In main.py:
# rmax has been set to 250au
# Added logic ~ If alpha=1e-3 set mmax=10, otherwise mmax=1e8

startingDir=100
n=150

# Stationary (100-111)
radialRes=200
velocity=0
for alpha in 3 4; do
  for amplitude in 3 10; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done

# Moving (112-123)
radialRes=240
position=90
for alpha in 3 4; do
  for amplitude in 3 10; do
    for velocity in 0.1 0.3 1; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done

# Miscellaneous prime scenarios
# Backwards
position=120
alpha=3
amplitude=10
velocity=-1
sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
((startingDir+=1))

