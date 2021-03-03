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
a3label = r"$\alpha = 10^{-3}$"
a4label = r"$\alpha = 10^{-4}$"
a3m = 'd'
a4m = '*'

# Center, ring width and fractional width data for moving bump
floor = 0
velocity = [1, 10, 30, 100]

# For SD > 10^-3
centerList310 = [108.1, 92, 51.2]
widthList310 = [12, 31.9, 64.3]
fracList310 = [0.11, 0.35, 1.26]
widthList410 = [4.2, 6.3, 13.1]
centerList410 = [114.1, 113.1, 109.2]
fracList410 = [0.04, 0.06, 0.13]

# New threshold rules
centerList310 = [118.4279, 109.204, 93.02, 53.6095]
widthList310 = [4.368, 14.08, 33.93, 69.08]
fracList310 = [0.03688, 0.12893, 0.36476, 1.2885]
widthList410 = [4.288, 6.374, 8.42, 16.2437]
centerList410 = [116.2636, 115.22, 114.197, 110.28587]
fracList410 = [0.03688, 0.055, 0.0737, 0.147288]

# Begin plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(12, 3.8))



ax1.set_ylabel("Ring width [au]", fontsize=fontsize)
ax1.plot(velocity, widthList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax1.plot(velocity, widthList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
ax1.legend(loc='upper left')
# ax1.set_xlim(0, 100)

ax2.set_ylabel("Ring center [au]", fontsize=fontsize)
ax2.loglog(velocity, centerList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax2.loglog(velocity, centerList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
# ax2.set_xlim(0, 100)
ax2.set_xlabel("Bump velocity $f$ as \% of nominal", fontsize=fontsize)

ax3.set_ylabel("Fractional width", fontsize=fontsize)
ax3.loglog(velocity, fracList310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label)
ax3.loglog(velocity, fracList410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label)
# ax3.set_xlim(0, 100)

# ax1.set_yscale('symlog')
# ax1.set_xscale('symlog')
# ax2.set_yscale('symlog')
# ax2.set_xscale('symlog')
# ax3.set_yscale('symlog')
# ax3.set_xscale('symlog')


# Formatting
fig.tight_layout()
filename = localDir + '/figplots/ringstatsloglog'
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(filename + '.png', dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
# plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()
