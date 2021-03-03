#!/bin/bash -l


#for dir in 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167; do
# 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201
for dir in 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217; do
  #python plotDistribution.py -z $dir -s 0
  #python plotSim.py -z $dir -s 1
  python plotd2g.py -z $dir
  #python plotMass.py -z $dir -s 0
done
#
#for dir in 138; do
#  python plotTime.py -z $dir
