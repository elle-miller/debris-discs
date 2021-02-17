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

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14


def main(args):

    # Overlay seaborn's styling with personal adjustments
    plt.style.use('seaborn-paper')
    plt.style.use('tex')
    plt.rcParams["figure.figsize"] = width_inches, height_inches

    # Read all data in the directory
    z = args.z
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
    center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
    textstr = getText(PlanMass[-1], center, width, frac)
    titlestr = getTitle(z, w)

    # Plot the surface density of dust and gas vs the distance from the star
    fig, ax = plt.subplots()
    ax.loglog(R[-1, ...], SigmaGas[-1, ...], label="Gas")
    ax.loglog(R[-1, ...], SigmaDustTot[-1, ...], label="Dust")
    ax.loglog(R[-1, ...], SigmaPlan[-1, ...], label="Plan")
    ax.loglog(R[-1, ...], d2g[-1, ...], ls='--', label="d2g")
    ax.set_ylim(1.e-6, 1.e4)
    ax.set_xlim(12, 200)
    ax.set_xlabel("Distance from star [au]")
    ax.set_ylabel("Surface density [g/cm²]")
    ax.legend(loc='upper right', frameon=True)
    ax.set_title(titlestr, fontdict={'fontsize': fontsize})
    ax.text(0.04, 0.85, textstr, transform=ax.transAxes)
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(("%.0f"))
    for axis in [ax.xaxis]:
        axis.set_major_formatter(ScalarFormatter())
        ax.set_xticks([30, 60, 90, 120])

    # Saving figure
    filename = outputDir + 'sdr/r' + str(z)
    fig.tight_layout()
    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
    if args.show:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)

