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
dirs = [259, 238, 208]
dirs = [277, 285, 284, 262, 208, 261]

labels = [r"$n = 0.03: \Delta r = 66$" + " au, " + r"$\Delta r/r = 1.3$", r"$n = 0.01: \Delta r = 55$" + " au, " + r"$\Delta r/r = 1.2$",
          r"${\rm step}\ \ \ \ \ \ : \Delta r$ = 45" + " au, " + r"$\Delta r/r = 1.1$"]
labels = [r"$n = 0.03$", r"$n = 0.01$", r"${\rm step}$"]

colors = ["C3", "C0", "C1"]
ls = ["-", "--", "-."]
ls = ["-", "-", "-", "-", "-", "-"]

f_vals = np.array([0.0, 0.01, 0.03])
n = len(f_vals)
n = 6
hex = [0xb7e6a5, 0x7ccba2, 0x46aea0, 0x089099, 0x00718b]
rgb = np.array([(252, 225, 164), (250, 191, 123), (240, 143, 110), (224, 92, 92), (209, 41, 89), (171, 24, 134)]) / 255
rgb = np.array([(252, 225, 164), (250, 191, 123), (240, 143, 110), (224, 92, 92), (209, 41, 89)]) / 255
rgb = np.array([(183, 230, 165), (124, 203, 162), (70, 174, 160), (8, 144, 153), (0, 113, 139), (4, 82, 117)]) / 255
# rgb = np.array([(158, 201, 226), (60, 147, 194), (13, 74, 112)]) / 255  # Mono-hue blue
# colors = np.array([rgb[5], rgb[3], rgb[1]])
colors = np.array(["midnightblue", "mediumseagreen", "lightsteelblue"])
colors = np.array(["tomato", "coral", "lemonchiffon"])
cmap = LinearSegmentedColormap.from_list("my_list", rgb, N=6)

i = 0
for z in dirs:

    # Read all data in the directory
    print("Sim #", z)
    w.datadir = getDataDir(z)
    R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    ax.loglog(R[-1, ...], SigmaPlan[-1, ...], ls=ls[i], color=cmap(i))
    i += 1


c1 = np.arange(1., n + 1)
norm = mplc.BoundaryNorm(np.arange(len(c1)+1)+0.5, len(c1))
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cbar = plt.colorbar(sm, ticks=c1, pad=.02)
#cbar.set_ticklabels([r"$\to 0$", r"$0.01$", r"$0.03$"])
cbar.set_ticklabels([r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$0.3$"])

# cbar.set_label(r"$n$", rotation=0, x=-0.5)
#fig.text(0.8, 0.05, '$n$', ha='center', va='center', rotation=0, fontsize=fontsize)
fig.text(0.8, 0.05, r'$\zeta$', ha='center', va='center', rotation=0, fontsize=fontsize)
cbar.ax.tick_params(which='both', size=0)

ax.set_xlabel("Distance from star [au]")
ax.set_ylabel("Surface density [g/cm²]")
filename = localDir + '/figplots/planratecmap'
filename = localDir + '/figplots/zeta'
# if args.title:
#     ax.set_title(titlestr, fontdict={'fontsize': fontsize})
#     filename += '_titled'
# if args.text:
#     ax.text(0.04, 0.85, textstr, transform=ax.transAxes)
# ax.xaxis.set_minor_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(("%.0f"))
# ax.yaxis.set_minor_formatter(("%.1f"))
# ax.yaxis.set_major_formatter(("%.1f"))
#ax.set_yticks([0.2, 0.1, 0.8, 0.6, 0.4])


ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(("%.0f"))
ax.set_xticks([30, 60])

# ax.yaxis.set_major_formatter(("%.2f"))
# ax.yaxis.set_major_formatter(("%.2f"))
ax.tick_params(axis='y', which='both', labelleft='on')
ax.get_yaxis().get_major_formatter().labelOnlyBase = False
# ax.set_yticks([0.1, 0.2])

# for axis in [ax.yaxis]:
#     axis.set_major_formatter(("%.2f"))
#     ax.set_yticks([0.2, 0.1, 0.04])
ax.set_ylim(0.02, 0.3)
ax.set_ylim(1e-3, 0.2)
ax.yaxis.set_minor_formatter(("%.2f"))
xmin = 17
xmax = 89
ax.set_xlim(xmin, xmax)
# maxPlan = np.max(SigmaPlan[-1])
# ax.hlines(maxPlan, xmin, xmax, 'k', label="Max")
# ax.hlines(0.75*maxPlan, xmin, xmax, 'r', label="0.75 Max")
# ax.hlines(0.5*maxPlan, xmin, xmax, 'g', label="0.50 Max")
# ax.hlines(0.25*maxPlan, xmin, xmax, 'm', label="0.25 Max")
# ax.hlines(0.1 * maxPlan, xmin, xmax, 'c', label="0.1 Max")


# Saving figure
#fig.tight_layout()
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
#plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])

plt.show()



