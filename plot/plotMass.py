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
slurmDir = '/mnt/beegfs/bachelor/groups/henning/users/miller/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14


def main(args):

    # Overlay seaborn's styling with personal adjustments
    plt.style.use('seaborn-paper')
    #plt.style.use('tex')
    plt.rcParams["figure.figsize"] = width_inches, height_inches
    fig, ax = plt.subplots()

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
    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
    SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    DustMass = (np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :])
    DustDiskMassEarth = np.sum(DustMass, axis=1) / M_earth
    PlanDiskMassEarth = w.read.sequence('planetesimals.M') / M_earth
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    rho_gas = w.read.sequence('gas.rho')
    rho_dust = w.read.sequence('dust.rho').sum(-1)
    d2g_mid = rho_dust / rho_gas
    DustVRad = w.read.sequence('dust.v.rad')
    OmegaK = w.read.sequence('grid.OmegaK')
    St = w.read.sequence('dust.St')

    # Text plotting
    RingDustTot = SigmaDustTot.copy()
    RingDiskMass = np.zeros(Nt)
    if z <= 201:
        [alpha0, amplitude, position] = getJobParams(z)
    else:
        position = 90
    iguess = (np.abs(R[0] - position)).argmin()

    # Obtain initial external dust
    it = 0
    initialRightDustTot = SigmaDustTot[it].copy()
    for i in range(Nr):
        if i <= iguess:
            initialRightDustTot[i] = 0
    initialRightDustMass = np.sum(np.pi * (rInt[it, 1:] ** 2. - rInt[it, :-1] ** 2.) * initialRightDustTot[:]) / M_earth
    d2g_mid_at_peak = np.zeros(Nt)
    dist = np.zeros(Nt)
    r_peak = np.zeros(Nt)
    M_flux = np.zeros(Nt)
    M_plan = np.zeros(Nt)
    SigmaDustTotPeak = np.zeros(Nt)
    imax = 141
    for i in range(Nt):
        igap = np.argmin(SigmaGas[i, 0:iguess + 10])
        iguess = igap
        ipeak = np.argmax(SigmaDustTot[i, igap:igap + 35]) + igap
        SigmaDustTotPeak[i] = SigmaDustTot[i, ipeak]
        dist[i] = R[i, ipeak] - R[i, igap]
        Rstart = R[i, ipeak] + 0.5 * dist[i]
        Rend = R[i, ipeak] + 0.5 - dist[i]
        istart = int(np.argmin(R[-1]-Rstart))
        iend = int(np.argmin(R[-1]-Rend))
        center, width, frac, istart2, iend2 = getRingStats(SigmaPlan[i], w)
        r_peak[i] = R[i, ipeak]
        d2g_mid_at_peak[i] = d2g_mid[i, ipeak]
        if ipeak != igap:
            RingDiskMass[i] = getRingMass(RingDustTot[i], istart, iend, w)


        # Calculate plan formation at each epoch
       #switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid_at_peak[i])) / 0.03))

        # Calculate inward planetesimal mass flux at each epoch
       #iflux = ipeak + (ipeak-igap)
       #flux = np.zeros(imax)
       #ret = np.zeros(imax)
       #for im in range(imax):
       #    flux[im] = SigmaDust[i, iflux, im] * DustVRad[i, iflux, im]
       #    ret[im] = 0.1 * SigmaDust[i, ipeak, im] * St[i, ipeak, im] * OmegaK[i, ipeak] * switch
       #total_flux = flux.sum()
       #dr = rInt[i, ipeak+1] - rInt[i, ipeak]
        #_plan[i] = 2 * c.pi * R[i, ipeak] * c.au * dr * ret.sum() / M_earth * c.year
       #M_flux[i] = -2 * c.pi * R[i, iflux] * c.au * np.sum(flux) / M_earth * c.year

    #ean_flux = np.mean(M_flux[-18:])
    #rint("Mean M_flux = ", mean_flux)
    textstr = getText(PlanDiskMassEarth[-1], center, width, frac, initialExtDust=initialRightDustMass)
    titlestr = getTitle(z, w)
    # ax.loglog(t, d2g_mid_at_peak, label="mid d2g peak")
    # ax.loglog(t, SigmaDustTotPeak, label="SigmaDust @ peak")
    # ax.loglog(t, M_flux, label="Mflux")
    # ax.hlines(mean_flux, 4e6, t[-1], label="mean ")

    ax.loglog(t, GasDiskMassEarth, label="Gas", color="C0")
    ax.loglog(t, DustDiskMassEarth, label="Dust", color="C1")
    ax.loglog(t, RingDiskMass, ls='--', label="Ring Dust", color="C4")
    ax.loglog(t, PlanDiskMassEarth, label="Planetesimals", color="C2")
    #ax.legend(loc='upper right')

    filename = outputDir + 'mass/m' + str(z)

    if args.title:
        ax.set_title(titlestr, fontdict={'fontsize': fontsize})
        filename += '_untitled'
    if args.text:
        ax.text(0.04, 0.85, textstr, transform=ax.transAxes)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Mass [M$_\oplus$]")
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:gray'
    ax2.set_ylabel('Midplane d2g ratio at peak', color=color, rotation=-90)  # we already handled the x-label with ax1
    ax2.loglog(t, d2g_mid_at_peak, '-.', linewidth=0.5, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    # Saving figure
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
    parser.add_argument('-t', action="store", dest="title", type=int, default=1, help="Show plot")
    parser.add_argument('-x', action="store", dest="text", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)
