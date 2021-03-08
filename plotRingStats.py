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
# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'
outputDir = localDir + '/figplots/'

width_inches = 6  # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
msize = 10

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
plt.style.use('tex')
# plt.rcParams["figure.figsize"] = width_inches, height_inches
a3col = 'darkmagenta'
a4col = 'orchid'
a3label = r"$\alpha = 10^{-3}$"
a4label = r"$\alpha = 10^{-4}$"
a3m = 'd'
a4m = '*'

# Center, ring width and fractional width data for moving bump
floor = 0
velocity = [10, 30, 100]

# Ring start and end information for moving bump
# 10 MYR
# s310 = [103.115, 76.764, 19.248]
# e310 = [115.182, 108.981, 87.343]
# s410 = [113.077, 111.010, 103.115]
# e410 = [117.326, 117.326, 117.326]
# timeLeft = (10 - 7.3564225445964215)  # Myr
# B = 14193.237814663926 # au/Myr
# distLeft = 1e-3*B*timeLeft
# s310[2] -= distLeft

dirs3 = np.array([206, 207, 208])
dirs4 = np.array([214, 215, 216])
s310 = np.zeros_like(dirs3, dtype=float)
e310 = np.zeros_like(dirs3, dtype=float)
w310 = np.zeros_like(dirs3, dtype=float)
c310 = np.zeros_like(dirs3, dtype=float)
f310 = np.zeros_like(dirs3, dtype=float)

s410 = np.zeros_like(dirs4, dtype=float)
e410 = np.zeros_like(dirs4, dtype=float)
w410 = np.zeros_like(dirs3, dtype=float)
c410 = np.zeros_like(dirs3, dtype=float)
f410 = np.zeros_like(dirs3, dtype=float)
for d in [dirs3, dirs4]:
    for i in range(len(d)):
        w.datadir = getDataDir(d[i])
        SigmaPlan = w.read.sequence('planetesimals.Sigma')
        R = w.read.sequence('grid.r') / c.au
        center, width, frac, istart, iend = getRingStats(SigmaPlan[-2], w)
        if d[i] == 208:
            center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
        if d is dirs3:
            s310[i] = float(R[-1, istart])
            e310[i] = float(R[-1, iend])
            c310[i] = center
            w310[i] = width
            f310[i] = frac
        else:
            s410[i] = R[-1, istart]
            e410[i] = R[-1, iend]
            c410[i] = center
            w410[i] = width
            f410[i] = frac

# # 7.26 MYR
# s310 = [117.326, 106.990, 87.343, 19.248]
# e310 = [119.510, 115.182, 108.981, 84.180]
# s410 = [115.182, 113.077, 111.010, 105.034]
# e410 = [117.326, 117.326, 117.326, 117.326]
#  #1e-4
# s310 = [117.326, 106.990, 85.747, 19.248]
# e310 = [119.510, 115.182, 108.981, 84.180]
# s410 = [115.182, 113.077, 111.010, 105.034]
# e410 = [117.326, 117.326, 117.326, 117.326]

# for i in range(3):
# c310 = np.array([0.5*(a + b) for a, b in zip(s310, e310)], dtype=float)
# w310 = np.array([b-a for a, b in zip(s310, e310)], dtype=float)
# f310 = np.array([a/b for a, b in zip(w310, c310)], dtype=float)
# c410 = np.array([0.5*(a + b) for a, b in zip(s410, e410)], dtype=float)
# w410 = np.array([b-a for a, b in zip(s410, e410)], dtype=float)
# f410 = np.array([a/b for a, b in zip(w410, c410)], dtype=float)

# print(w310)
# print(w410)

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
endtime = 7.3564225445964215

plotPositions3 = False
if plotPositions3:

    fig, ax = plt.subplots()
    ax.plot(velocity, s310, markersize=msize, marker='*', ls='--', color="C0", label=r"$\alpha_0$ = 1e-3")
    ax.plot(velocity, e310, markersize=msize, marker='*', ls='--', color="C0")
    ax.plot(velocity, s410, markersize=msize, marker='.', ls='-.', color="C3", label=r"$\alpha_0$ = 1e-4")
    ax.plot(velocity, e410, markersize=msize, marker='.', ls='-.', color="C3")
    ax.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)
    ax.set_ylabel("Final ring location [au]", fontsize=fontsize)
    ax.legend(fontsize=fontsize, loc='lower left')
    ax.set_ylim(0, 140)
    ax.set_title("A = 10 bump at 90 au evolved for 10 Myr", fontsize=fontsize)
    fig.tight_layout()
    plt.savefig(outputDir + 'ringpositions.png', dpi=300)
    plt.show()

plotPositionTime = True
if plotPositionTime:

    alpha3 = True
    if alpha3:
        dirs = [208, 207, 206]
        filename = 'ringevolution4.png'
        ymin = 15
    else:
        dirs = [216, 215, 214]
        filename = 'ringevolution3.png'
        ymin = 70
    f_vals = np.array([100, 30, 10])
    f = [r"$f = 10\%$", r"$f = 30\%$", r"$f = 100\%$"]
    colors = ["C0", "C1", "C3"]
    facecolors = ["powderblue", "navajowhite", "lightcoral"]
    fig, ax = plt.subplots()

    n = len(f_vals)
    hex = [0xFCE1A4, 0xFABF7B, 0xF08F6E, 0xE05C5C, 0xd12959, 0xab1886]
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
    cbar.ax.set_ylabel(r"Bump velocity $f$")
    cbar.ax.tick_params(which='both', size=0)

    ax.plot(t, np.zeros_like(rgap), ls='--', color='lightgray', label='Bump position', linewidth=1.5)
    ax.legend(fontsize=fontsize, loc='lower left')

    ax.set_xlim(3, 7.35)
    ax.set_xlabel("Time [Myr]", fontsize=fontsize)
    ax.set_ylabel("Belt position [au]", fontsize=fontsize)
    ax.set_ylim(ymin, 125)
    # ax.set_title(titlestr, fontdict={'fontsize': fontsize})
    plt.savefig(outputDir + filename, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.show()
    exit(0)


yeet = True
if yeet:
    # Begin plotting
    fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))

    plt.subplots_adjust(hspace=0.1)
    ax1.set_ylabel("Belt position [au]", fontsize=fontsize)
    size=0.9
    velocity = [10, 30, 100]
    ax1.plot(velocity, s310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label, linewidth=size)
    ax1.plot(velocity, e310, markersize=msize, marker=a3m, ls='--', color=a3col, linewidth=size)
    ax1.plot(velocity, s410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label, linewidth=size)
    ax1.plot(velocity, e410, markersize=msize, marker=a4m, ls='--', color=a4col, linewidth=size)
    ax1.tick_params(axis='both', which='both', top='on', bottom='on', right='on')
    ax1.legend(loc='lower left')

    ax3.set_xlabel("Bump velocity $f$ as \% of nominal", fontsize=fontsize)
    ax3.set_ylabel("Fractional width", fontsize=fontsize)
    ax3.plot(velocity, f310, markersize=msize, marker=a3m, ls='--', color=a3col, label=a3label, linewidth=size)
    ax3.plot(velocity, f410, markersize=msize, marker=a4m, ls='--', color=a4col, label=a4label, linewidth=size)
    ax3.tick_params(axis='both', which='both', top='on', bottom='on', right='on')


    filename = localDir + '/figplots/ringstats10Myr'
    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(filename + '.png', dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
    # plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
    plt.show()
