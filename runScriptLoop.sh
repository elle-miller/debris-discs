#!/bin/bash -l


for dir in 70 71 72 73 74 75 76 77 78 79 80 81; do
  python getplanmasses.py -a 1 -z $dir
  #python plotSim.py -r 1 -z $dir
done
