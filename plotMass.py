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
from plottingFunctions import getDataDir, getTitle, getText, getRingStats, getRingMass
from dustpy import hdf5writer as w
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
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
    rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
    SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    DustMass = (np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :])
    DustDiskMassEarth = np.sum(DustMass, axis=1) / M_earth
    PlanDiskMassEarth = w.read.sequence('planetesimals.M') / M_earth
    SigmaPlan = w.read.sequence('planetesimals.Sigma')

    # Text plotting
    RingDustTot = SigmaDustTot.copy()
    RingDiskMass = np.zeros(Nt)
    for i in range(Nt):
        center, width, frac, istart, iend = getRingStats(SigmaPlan[i], w)
        RingDiskMass[i] = getRingMass(RingDustTot[i], istart, iend, w)
    textstr = getText(PlanDiskMassEarth[-1], center, width, frac)
    titlestr = getTitle(z, w)

    fig, ax = plt.subplots()
    ax.loglog(t, GasDiskMassEarth, label="Gas", color="C0")
    ax.loglog(t, DustDiskMassEarth, label="Dust", color="C1")
    ax.loglog(t, RingDiskMass, ls='--', label="Ring Dust", color="C4")
    ax.loglog(t, PlanDiskMassEarth, label="Planetesimals", color="C2")
    ax.set_xlim(1e4, 1e7)
    ax.set_ylim(1e-1, 1e7)
    ax.legend(loc='upper right')
    ax.lineTime = ax.axvline(t[0], color="C7", zorder=-1, lw=1)
    ax.set_title(titlestr, fontdict={'fontsize': fontsize})
    ax.text(0.04, 0.85, textstr, transform=ax.transAxes)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Mass [M$_\oplus$]")

    # Saving figure
    filename = outputDir + 'mass/m' + str(z)
    fig.tight_layout()
    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(p["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    plt.savefig(e["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    if args.show:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)
