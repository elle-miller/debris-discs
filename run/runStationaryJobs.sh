#!/bin/bash -l

# Bash script to run the mpi script multiple times for varying parameters
# Date: Final rerun of stationary scripts with higher Nt

startingDir=1
n=150
radialRes=200
for alpha in 3 4; do
  for amplitude in 3 10; do
    for position in 30 60 90; do
      sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
      ((startingDir+=1))
    done
  done
done


#alpha=3
#amplitude=10
#position=90
#velocity=0
#for v in 0 1; do
#  sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#  ((startingDir+=1))
#done




## 190 double
#startingDir=240
#n=31
#radialRes=200
#alpha=3
#amplitude=30
#position=30
#velocity=0
#sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes
#
## 210 double
#startingDir=241
#n=31
#radialRes=200
#alpha=4
#amplitude=3
#position=90
#velocity=0.1
#sbatch job.sh $startingDir $alpha $amplitude $velocity $position $n $radialRes

