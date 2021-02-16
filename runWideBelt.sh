#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 16 - create a wide belt based on 208

startingDir=232
startyear=5
velocity=1
alpha=3
amplitude=10

# Case 1
position=120
timemove=0
sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear
((startingDir+=1))

# Case 2
position=90
timemove=1000000
sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear
((startingDir+=1))

# Case 3
position=120
timemove=1000000
sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear


