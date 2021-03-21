#!/usr/bin/python3
# -*- coding: utf-8 -*-

from dustpy import plot
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import numpy as np
import argparse
import matplotlib.pyplot as plt
from movieBump import movieBump
from dustpy.plot import panel
from os import path, getcwd
from matplotlib.ticker import ScalarFormatter
from plottingFunctions import *
from dustpy import hdf5writer as w
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mplc

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
outputDir = localDir + '/simplots/'

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches

fig, ax = plt.subplots()
dirs = [208, 272]

labels = ["90 au", "85 au"]

colors = ["C3", "C0", "C1"]
ls = ["-", "--", "-."]
ls = ["-", "-", "-"]

f_vals = np.array([0.0, 0.01, 0.03])
n = len(f_vals)
hex = [0xb7e6a5, 0x7ccba2, 0x46aea0, 0x089099, 0x00718b]
rgb = np.array([(252, 225, 164), (250, 191, 123), (240, 143, 110), (224, 92, 92), (209, 41, 89), (171, 24, 134)]) / 255
# rgb = np.array([(183, 230, 165), (124, 203, 162), (70, 174, 160), (8, 144, 153), (0, 113, 139), (4, 82, 117)]) / 255
# rgb = np.array([(158, 201, 226), (60, 147, 194), (13, 74, 112)]) / 255  # Mono-hue blue
colors = np.array([rgb[5], rgb[3], rgb[1]])
# colors = np.array(["midnightblue", "royalblue", "lightsteelblue"])
cmap = LinearSegmentedColormap.from_list("my_cmap", colors, N=n)

i = 0
for z in dirs:

    # Read all data in the directory
    print("Sim #", z)
    w.datadir = getDataDir(z)
    R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
    t = w.read.sequence('t')
    print(t[-1] / c.year * 1e-6)
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    ax.loglog(R[-1, ...], SigmaPlan[-1, ...], ls=ls[i], color=colors[i], label=labels[i])
    i += 1


ax.set_xlabel("Distance from star [au]")
ax.set_ylabel("Surface density [g/cm²]")
filename = localDir + '/figplots/coag'


ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(("%.0f"))
ax.set_xticks([30, 60])

# ax.yaxis.set_major_formatter(("%.2f"))
# ax.yaxis.set_major_formatter(("%.2f"))
ax.tick_params(axis='y', which='both', labelleft='on')
ax.get_yaxis().get_major_formatter().labelOnlyBase = False
# ax.set_yticks([0.1, 0.2])
ax.legend()

ax.set_ylim(6e-3, 1)
ax.yaxis.set_minor_formatter(("%.2f"))
xmin = 10
xmax = 89
ax.set_xlim(xmin, xmax)


e = getEPS(filename)
p = getPNG(filename)
plt.savefig(outputDir + "diffinitr.png", dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
# plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])

plt.show()



