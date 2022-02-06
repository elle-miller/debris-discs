#!/bin/bash -l

startingDir=31
alpha=3
amplitude=10
velocity=0
position=90
radialRes=240
n=31

sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
