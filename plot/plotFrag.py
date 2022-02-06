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
msize = 8

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches
a3col = 'lightcoral'
a4col = 'cornflowerblue'

# Center, ring width and fractional width data for moving bump
velocity = [1, 3, 10]
mass30 = [98, 213, 238]
mass60 = [1, 110, 135]
mass90 = [0, 55, 72]

# Begin plotting
fig, ax1 = plt.subplots(1, 1)
ax1.set_ylabel("Ring width [au]", fontsize=fontsize)
ax1.plot(velocity, mass30, markersize=msize, marker='o', ls='--', label=r"r$_{\rm p} = 30$ au")
ax1.plot(velocity, mass60, markersize=msize, marker='^', ls='--', label=r"r$_{\rm p} = 60$ au")
ax1.plot(velocity, mass90, markersize=msize, marker='s', ls='--', label=r"r$_{\rm p} = 90$ au")
# ax1.semilogx(velocity, mass30, markersize=msize, marker='o', ls='--', label=r"r$_{\rm p} = 30$ au")
# ax1.semilogx(velocity, mass60, markersize=msize, marker='^', ls='--', label=r"r$_{\rm p} = 60$ au")
# ax1.semilogx(velocity, mass90, markersize=msize, marker='s', ls='--', label=r"r$_{\rm p} = 90$ au")
ax1.legend()
ax1.set_xlabel("Fragmentation velocity [m/s]", fontsize=fontsize)
ax1.set_ylabel("Planetesimal mass [M$_\oplus$]", fontsize=fontsize)

# Formatting
fig.tight_layout()
filename = localDir + '/figplots/vfrag'
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()


# Alternative
rp = [30, 60, 90]
v1 = [98, 1, 0]
v3 = [213, 110, 55]
v10 = [238, 135, 72]

fig, ax1 = plt.subplots(1, 1)
ax1.set_ylabel("Ring width [au]", fontsize=fontsize)
# ax1.plot(rp, v1, markersize=msize, marker='o', ls='--', label=r"v$_{\rm f} = 1$ m/s")
# ax1.plot(rp, v3, markersize=msize, marker='^', ls='--', label=r"v$_{\rm f} = 3$ m/s")
# ax1.plot(rp, v10, markersize=msize, marker='s', ls='--', label=r"v$_{\rm f} = 10$ m/s")


ax1.semilogx(rp, v10, markersize=msize, marker='s', ls='--', label=r"v$_{\rm f} = 10$ m/s")
ax1.semilogx(rp, v3, markersize=msize, marker='^', ls='--', label=r"v$_{\rm f} = 3$ m/s")
ax1.semilogx(rp, v1, markersize=msize, marker='o', ls='--', label=r"v$_{\rm f} = 1$ m/s")
ax1.set_ylim(3, 250)
ax1.xaxis.set_major_formatter(ScalarFormatter())
ax1.xaxis.set_minor_formatter(ScalarFormatter())
# ax1.xaxis.set_major_formatter(("%.0f"))
ax1.legend()
ax1.set_xlabel(r"r$_{\rm p}$ [au]", fontsize=fontsize)
ax1.set_ylabel("Planetesimal mass [M$_\oplus$]", fontsize=fontsize)

# Formatting
fig.tight_layout()
filename = localDir + '/figplots/vfragpos'
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()
