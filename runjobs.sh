#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: October 20
alpha=4
startingDir=360


# Two areas where alpha is a big problem
# Stationary bump: 1e-4, A=3, pos=30
A1=3
v1=0
p1=30

# Moving bump: 1e-4, A=10, v=0.1, pos=90
A2=10
v2=0.1
p2=90

n=3

# Stationary bump
sbatch job.sh $startingDir $alpha $A1 $v1 $p1 $n "$0"
((startingDir+=1))
# Moving bump
sbatch job.sh $startingDir $alpha $A2 $v2 $p2 $n "$1"



## Loop through each possible combination, calling job.sh
#Nr=1000
#for dustSolver in "IMPL_DIRECT" "EXPL_EULER"  ; do
#  # Stationary bump
#  sbatch job.sh $startingDir $alpha $A1 $v1 $p1 $Nr $dustSolver
#  ((startingDir+=1))
#  # Moving bump
#  sbatch job.sh $startingDir $alpha $A2 $v2 $p2 $Nr $dustSolver
#  ((startingDir+=1))
#done
