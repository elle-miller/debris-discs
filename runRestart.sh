#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Feb 23

# 190
startingDir=190
amplitude=30
position=30
velocity=0
sbatch jobRestart.sh $startingDir $amplitude $velocity $position

# 210
startingDir=210
amplitude=3
position=30
velocity=0.1
sbatch jobRestart.sh $startingDir $amplitude $velocity $position
