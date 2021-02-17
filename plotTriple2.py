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

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14

# Read all data in the directory
fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(width_inches, 3.5*height_inches))
fontsize = 14

# Overlay seaborn's styling with personal adjustments
plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches

plot3 = True
if plot3:
    firstendring = 115
    dir = [206, 207, 208]
    filename = localDir + '/figplots/movingcompare3'
else:
    firstendring = 117
    dir = [214, 215, 216]
    filename = localDir + '/figplots/movingcompare4'
i = 0
for z in dir:

    # Read all data in the directory
    print("Sim #", z)
    w.datadir = getDataDir(z)
    outputDir = localDir + '/simplots/'

    # Data
    R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    d2g = w.read.sequence('dust.eps')
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
    PlanMass = w.read.sequence('planetesimals.M') / M_earth

    # Text plotting
    [alpha0, amplitude, position] = getJobParams(z)
    center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
    textstr = getText(PlanMass[-1], center, width, frac)
    titlestr = "f = {p}%".format(p=position)

    # Plot the surface density of dust and gas vs the distance from the star
    ax[i].loglog(R[-1, ...], SigmaGas[-1, ...], label="Gas")
    ax[i].loglog(R[-1, ...], SigmaDustTot[-1, ...], label="Dust")
    ax[i].loglog(R[-1, ...], SigmaPlan[-1, ...], label="Plan")
    ax[i].loglog(R[-1, ...], d2g[-1, ...], ls='--', label="d2g")
    ax[i].set_ylim(1.e-6, 1.e4)
    ax[i].set_xlim(12, 200)
    ax[i].set_title(titlestr, fontdict={'fontsize': fontsize, 'font': 'serif'})
    ax[i].text(0.04, 0.88, textstr, transform=ax[i].transAxes)
    ax[i].xaxis.set_major_formatter(ScalarFormatter())
    ax[i].vlines(firstendring, 1.e-6, 1e4, 'r')
    i = i + 1

fig.text(0.5, 0.01, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
fig.text(0.02, 0.5, 'Surface density [g/cm²]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
fig.tight_layout()
e = getEPS(filename)
p = getPNG(filename)
plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
plt.show()