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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mplc
from matplotlib.cm import get_cmap
from matplotlib.legend_handler import HandlerLine2D
# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
slurmDir = '/mnt/beegfs/bachelor/groups/henning/users/miller/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'
outputDir = localDir + '/figplots/'

width_inches = 6  # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
msize = 8

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
# plt.style.use('tex.mplstyle')
plt.rcParams["figure.figsize"] = width_inches, height_inches
plt.rcParams["axes.spines.left"] = True  # display axis spines
plt.rcParams["axes.spines.bottom"] = True
plt.rcParams["axes.spines.top"] = True
plt.rcParams["axes.spines.right"] = True
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["font.size"] = 14
plt.rcParams["legend.fontsize"] = 14
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["legend.labelspacing"] = 0.2
# Read all data in the directory
outputDir = localDir + '/figplots/'

# Center, ring width and fractional width data for moving bump
floor = 0
velocity = [0, 10, 30, 100]

# constants
au = 1.495978707e11
k = 1.380649e-23
mu = 2.3
mp = 1.67262192369e-27
G = 6.6743e-11
M = 1.988409870698051e30  # kg
R = 2 * 1391400000  # stellar radius (m)
T = 5772
phi = 0.05
B = (1.5 * k * T)/(mu * mp) * (R * np.sqrt(phi)/(G * M))**0.5
B = B / au * c.year
alpha3 = 1e-3
alpha4 = 1e-4
endtime = 7.356422544596421
endtime = 7.57172149

alpha3 = True
if alpha3:
    dirs = [117, 116, 115]
    filename = 'ringevolution4_new'
    ymin = 15
else:
    dirs = [216, 215, 214]
    filename = 'ringevolution3_new'
    ymin = 70
f_vals = np.array([100, 30, 10])
f = [r"$f = 10\%$", r"$f = 30\%$", r"$f = 100\%$"]
fig, ax = plt.subplots()

n = len(f_vals)
rgb = np.array([(252, 225, 164), (250, 191, 123), (240, 143, 110), (224, 92, 92), (209, 41, 89), (171, 24, 134)]) / 255
colors = np.array([rgb[1], rgb[3], "firebrick"], dtype=object)
cmap = LinearSegmentedColormap.from_list("my_cmap", colors, N=n)

for i in range(len(dirs)):
    z = dirs[i]
    w.datadir = getDataDir(z)
    t = w.read.sequence('t') / c.year * 1e-6
    Nt = t.shape[0]
    R = w.read.sequence('grid.r') / c.au
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    titlestr = getTitle(z, w)
    rstart = np.zeros(Nt)
    rend = np.zeros(Nt)
    rgap = np.zeros(Nt)
    rpeak = np.zeros(Nt)

    iguess = (np.abs(R[0] - 90)).argmin()
    for it in range(Nt):
        igap = np.argmin(SigmaGas[it, 0:iguess + 1])
        iguess = igap
        ipeak = np.argmax(SigmaDustTot[it, igap:igap + 35]) + igap
        dist = ipeak - igap
        center, width, frac, istart, iend = getRingStats(SigmaPlan[it], w)
        rstart[it] = R[it, istart]
        rend[it] = R[it, iend]
        rgap[it] = R[it, igap]
        rpeak[it] = R[it, ipeak]

    ax.fill_between(t, rstart, rend, facecolor=cmap(i))
    ax.plot(t, rstart, markersize=msize, ls='-', color=cmap(i))
    ax.plot(t, rend, markersize=msize, ls='-', color=cmap(i))
    ax.plot(t, rgap,  ls='--', color=cmap(i), linewidth=1.5, zorder=10)

c1 = np.arange(1., n + 1)
norm = mplc.BoundaryNorm(np.arange(len(c1)+1)+0.5, len(c1))
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cbar = plt.colorbar(sm, ticks=c1, pad=.02)
cbar.set_ticklabels([r"$100\%$", r"$30\%$", r"$10\%$"])
cbar.ax.tick_params(which='both', size=0)
fig.text(0.8, 0.05, r'$f$', ha='center', va='center', rotation=0, fontsize=fontsize)

ax.plot(t, np.zeros_like(rgap), ls='--', color='lightgray', label='Gap position', linewidth=1.5)
ax.legend(fontsize=fontsize, loc='lower left')

ax.set_xlim(3, endtime)
ax.set_xlabel("Time [Myr]", fontsize=fontsize)
ax.set_ylabel("Belt position [au]", fontsize=fontsize)
ax.set_ylim(ymin, 120)
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(filename + '.png', dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
plt.savefig(filename + '.eps', dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.savefig(filename + '.pdf', dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()
exit(0)
