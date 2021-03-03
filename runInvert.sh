#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Sharp formation, tanh commented out
startingDir=266
n=31
radialRes=200
position=90
amplitude=10
velocity=1
alpha=3
invert=1
sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes $invert