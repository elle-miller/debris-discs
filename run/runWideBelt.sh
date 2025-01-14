#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: May 4 be with you - last experiments

startingDir=30
startyear=5
velocity=1
alpha=3
amplitude=10
Nr=240
n=31

# Case 1
# position=120
# timemove=0
# sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear $Nr $n
# ((startingDir+=1))

# Case 2
# position=90
# timemove=1000000
# sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear $Nr $n
# ((startingDir+=1))

# Case 3
# position=120
# timemove=1000000
# sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear $Nr $n

# Case 4
position=120
timemove=2000000
sbatch jobWideBelt.sh $startingDir $alpha $amplitude $velocity $position $timemove $startyear $Nr $n


