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
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
    Nr = R.shape[1]
    rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    d2g = w.read.sequence('dust.eps')
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
    PlanMass = w.read.sequence('planetesimals.M') / M_earth
    rho_gas = w.read.sequence('gas.rho')
    rho_dust = w.read.sequence('dust.rho').sum(-1)
    d2g_mid = rho_dust / rho_gas

    # Text plotting
    filename = outputDir + 'd2g/d2g' + str(z)
    center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
    # print("Width = ", width, ", center = ", center, ", Frac = ", frac)
    textstr = getText(PlanMass[-1], center, width, frac, justMass=True)
    titlestr = getTitle(z, w)

    fig, ax = plt.subplots()
    d2g_mid_at_peak = np.zeros(Nt)
    r_peak = np.zeros(Nt)
    # Starting guess
    if z <= 201:
        [alpha0, amplitude, position] = getJobParams(z)
    else:
        position = 90
    iguess = np.argmin(abs(R[-1] - position))
    for it in range(Nt):
        igap = np.argmin(SigmaGas[it, 0:iguess + 10])
        iguess = igap
        ipeak = np.argmax(SigmaDustTot[it, igap:igap + 35]) + igap
        r_peak[it] = R[it, ipeak]
        d2g_mid_at_peak[it] = d2g_mid[it, ipeak]
    #ax.loglog(t, r_peak)
    ax.loglog(t, d2g_mid_at_peak, 'r')
    # ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Time [yr]")
    ax.set_ylabel("Midplane dust-to-gas ratio at peak")

    if args.title:
        ax.set_title(titlestr, fontdict={'fontsize': fontsize})
        filename += '_titled'
    if args.text:
        ax.text(0.04, 0.85, textstr, transform=ax.transAxes)

    # Saving figure
    fig.tight_layout()
    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
    # plt.savefig(e["filename"], dpi=300, bbox_inches=e["bbox"], pad_inches=e["pad"])
    # if args.show:
    #     plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    parser.add_argument('-t', action="store", dest="title", type=int, default=1, help="Show plot")
    parser.add_argument('-x', action="store", dest="text", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)

