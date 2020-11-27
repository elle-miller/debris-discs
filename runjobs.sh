#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Nov 8 - running on Nov19
alpha=4
startingDir=2

# Two areas where alpha is a big problem
# Stationary bump: 1e-4, A=3, pos=30
A1=3
v1=0
p1=30

# Moving bump: 1e-4, A=10, v=0.1, pos=90
A2=10
v2=0.1
p2=90

# Planetesimal scripts with only 3 snapshots
n=3
planOff=0
planOn=1
# Stationary bump with formation
sbatch job.sh $startingDir $alpha $A1 $v1 $p1 $n $planOn
((startingDir+=1))

# Stationary bump with no formation
sbatch job.sh $startingDir $alpha $A1 $v1 $p1 $n $planOff
((startingDir+=1))

# Loop through two alphas, two amplitudes and three positions. 50 snapshots. PF ON
n=50
velocity=0
for alpha in 3 4; do
  for amplitude in 10 30; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $planOn
      ((startingDir+=1))
    done
  done
done
