#!/bin/bash -l

# Job with step planetesimal function
startingDir=239
n=31
radialRes=200
alpha=3
amplitude=10
position=90
velocity=1
sbatch jobFrag.sh "$startingDir" $alpha $amplitude $velocity $position $n $radialRes
