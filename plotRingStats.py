from dustpy import plot
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import numpy as np
import argparse
import matplotlib.pyplot as plt
from movieBump import movieBump
from dustpy.plot import panel
from os import path
from os import getcwd
from matplotlib.ticker import ScalarFormatter
from plottingFunctions import *

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6  # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
msize = 10

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
plt.style.use('tex')
# plt.rcParams["figure.figsize"] = width_inches, height_inches
a3col = 'lightcoral'
a4col = 'cornflowerblue'
a3label = r"$\alpha$ = 1e-3"
a4label = r"$\alpha$ = 1e-4"
a3m = 'd'
a4m = '*'

# Center, ring width and fractional width data for moving bump
floor = 0
velocity = [10, 30, 100]
centerList310 = [108.1, 92, 51.2]
widthList310 = [12, 31.9, 64.3]
fracList310 = [0.11, 0.35, 1.26]
widthList410 = [4.2, 6.3, 13.1]
centerList410 = [114.1, 113.1, 109.2]
fracList410 = [0.04, 0.06, 0.13]

# Begin plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(12, 3.8))
ax1.set_ylabel("Ring width [au]", fontsize=fontsize)
loc = 'upper left'
figname = 'width34'
ax1.plot(velocity, widthList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax1.plot(velocity, widthList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
ax1.legend(loc=loc)
#ax1.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)

ax2.set_ylabel("Ring center [au]", fontsize=fontsize)
loc = 'lower left'
figname = 'center34'
ax2.plot(velocity, centerList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax2.plot(velocity, centerList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
#ax2.legend(loc=loc)
ax2.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)

ax3.set_ylabel("Fractional width", fontsize=fontsize)
loc = 'upper left'
figname = 'frac34'
ax3.plot(velocity, fracList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax3.plot(velocity, fracList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
#ax3.legend(loc=loc)
#ax3.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)

# Formatting
fig.tight_layout()
filename = localDir + '/figplots/ringstats'
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()
