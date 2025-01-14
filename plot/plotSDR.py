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
# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
sdr_colors = [i for i in get_cmap('tab20').colors]

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

    # b = (R[-1, np.argmax(SigmaPlan[-1])])
    # a = (R[-1, np.argmax(SigmaPlan[30, 59])])

    # Text plotting
    filename = outputDir + 'sdr/r' + str(z)
    center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
    print("Width = ", width, ", center = ", center, ", Frac = ", frac)
    av_plan_sd = np.mean(SigmaPlan[-1, istart:iend])
    print("Average plan SD = ", av_plan_sd)
    textstr = getText(PlanMass[-1], center, width, frac)
    titlestr = getTitle(z, w)

    # Plot the surface density of dust and gas vs the distance from the star
    fig, ax = plt.subplots()
    ax.loglog(R[-1, ...], SigmaGas[-1, ...], color=sdr_colors[0], label="Gas")
    ax.loglog(R[-1, ...], SigmaDustTot[-1, ...], color=sdr_colors[2], label="Dust")
    ax.loglog(R[-1, ...], SigmaPlan[-1, ...], color=sdr_colors[4], label="Planetesimals",)
    # ax.loglog(R[-1, ...], d2g[-1, ...], color=sdr_colors[-2], ls='--', label="d2g")
    ax.loglog(R[-1, ...], d2g_mid[-1, ...], color=sdr_colors[6+8], ls='-.', label=r"$\rho_d/\rho_g$", linewidth=1.2)
    # ax.hlines(np.max(SigmaPlan[-1]), 1e-10, 1e4)
    #ax.hlines(8.5e-2, 1e-10, 1e4, label='8.5')
    #ax.hlines(av_plan_sd, 1e-10, 1e4, label='av')
    # ax.vlines(19.6, 1e-10, 1e4)
    # ax.vlines(31, 1e-10, 1e4)
    # ax.vlines(50, 1e-10, 1e4)
    ax.set_ylim(1e-5, 1e2)
    xmin = 10
    xmax = 200
    # xmax = 300
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Distance from star [au]")
    ax.set_ylabel("Surface density [g/cm²]")
    ax.legend(loc='center left')

    if args.title:
        ax.set_title(titlestr, fontdict={'fontsize': fontsize})
        filename += '_titled'
    if args.text:
        ax.text(0.04, 0.9, textstr, transform=ax.transAxes)
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(("%.0f"))
    for axis in [ax.xaxis]:
        axis.set_major_formatter(ScalarFormatter())
        ax.set_xticks([30, 60, 90, 120, 200])
    # Saving figure
    fig.tight_layout()
    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(p["filename"], dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])
    plt.savefig(filename+".pdf", dpi=300, bbox_inches=p["bbox"], pad_inches=p["pad"])

    # plt.savefig(e["filename"], dpi=300) #, bbox_inches=e["bbox"], pad_inches=e["pad"])
    if args.show:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    parser.add_argument('-t', action="store", dest="title", type=int, default=1, help="Show plot")
    parser.add_argument('-x', action="store", dest="text", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)

