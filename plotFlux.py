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
from matplotlib.cm import get_cmap

sdr_colors = [i for i in get_cmap('tab20').colors]


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
    fig, ax = plt.subplots()

    # Read all data in the directory
    z = args.z
    print("Sim #", z)
    w.datadir = getDataDir(z)
    outputDir = localDir + '/figplots/'

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
    M_plan = np.zeros([Nt, Nr])
    predicted_sd = np.zeros(Nt)
    centers = np.zeros(Nt)
    SigmaDustTotPeak = np.zeros(Nt)
    SigmaPlanPeak = np.zeros(Nt)
    imax = 141
    v_bump = -4757.6  # cm/s  10 * c.au / c.year * 1e-6  # cm/s

    for i in range(Nt):
        igap = np.argmin(SigmaGas[i, 0:iguess + 10])
        iguess = igap
        ipeak = np.argmax(SigmaDustTot[i, igap:igap + 35]) + igap
        SigmaDustTotPeak[i] = SigmaDustTot[i, ipeak]
        SigmaPlanPeak[i] = SigmaPlan[i, ipeak]
        dist[i] = R[i, ipeak] - R[i, igap]
        Rstart = R[i, ipeak] + 0.5 * dist[i]
        Rend = R[i, ipeak] + 0.5 - dist[i]
        istart = int(np.argmin(R[-1]-Rstart))
        iend = int(np.argmin(R[-1]-Rend))
        center, width, frac, istart2, iend2 = getRingStats(SigmaPlan[i], w)
        centers[i] = center
        r_peak[i] = R[i, ipeak]
        d2g_mid_at_peak[i] = d2g_mid[i, ipeak]
        if ipeak != igap:
            RingDiskMass[i] = getRingMass(RingDustTot[i], istart, iend, w)

        # Calculate plan formation at each epoch
        switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid_at_peak[i])) / 0.03))

        # Calculate inward planetesimal mass flux at each epoch
        iflux = (ipeak + 1) #(ipeak-igap))
        flux = np.zeros(imax)
        ret = np.zeros([imax, Nr])
        #print(DustVRad[-1, iflux, :])
        for im in range(imax):
            flux[im] = SigmaDust[i, iflux, im] * (DustVRad[i, iflux, im] - v_bump)
            ret[im, :] = (0.1 * SigmaDust[i, :, im] * St[i, ipeak, im] * OmegaK[i, :] * switch)
        #print(DustVRad[-1])

        total_flux = flux.sum()
        dr = rInt[i, ipeak+1] - rInt[i, ipeak]
        M_plan[i, :] = (2 * c.pi * R[i, :] * c.au * dr * np.sum(ret, axis=0) / M_earth * c.year)
        M_flux[i] = -2 * c.pi * R[i, iflux] * c.au * np.sum(flux) / M_earth * c.year
        predicted_sd[i] = M_flux[i] * M_earth / c.year / (2 * c.pi * centers[i] * c.au * v_bump)  # g s-1 cm-1 cm-1 s = g / cm2

    mean_flux = np.mean(M_flux[-100:])
    print("Mean M_flux = ", mean_flux)
    #textstr = getText(PlanDiskMassEarth[-1], center, width, frac, initialExtDust=initialRightDustMass)
    #titlestr = getTitle(z, w)
    #print(M_flux)

    lns1 = ax.semilogy(t * 1e-6, M_flux, color="C0", label=r"$\dot{M}_{\rm peb}$")
    lns2 = ax.semilogy(t * 1e-6, np.sum(M_plan, axis=1), color=sdr_colors[4], label=r"$\dot{M}_{\rm plan}$")
    #lns3 = ax.semilogy([1, 2], [1e-10, 1e-10], ls='--', color=sdr_colors[5], label=r"min($\Sigma_{\rm plan}$)")
    #ax.vlines(4.7, 1e-10, 1e-3, ls='--', color=sdr_colors[5])
    # ax.vlines(6, 1e-10, 1e-3, ls='--', color=sdr_colors[5])
    # ax.hlines(mean_flux, minyear, maxyear, label="mean ")
    ax.set_xlim(1, 8)
    ax.set_ylim(1e-9, 1e-2)
    filename = outputDir + "208mfluxplus1"
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Mass rate [M$_\oplus$/yr]")

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:gray'
    lns5 = ax2.plot(t * 1e-6, d2g_mid_at_peak, '-.', label=r"$\rho_d/\rho_g(r_{\rm peak})$", linewidth=1.5, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim(0.6, 1)

    lns = lns1 + lns2 +  lns5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='lower right', frameon=True, framealpha=0.8)

    e = getEPS(filename)
    p = getPNG(filename)
    plt.savefig(p["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    plt.savefig(e["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    if args.show:
        plt.show()

    exit(0)
    fig, ax = plt.subplots()
    ax.loglog(t, predicted_sd, label="predicted")
    ax.loglog(t, SigmaPlanPeak, label="real plan")
    ax.loglog(t, d2g_mid_at_peak, label=r"$\rho_d/\rho_g$")
    ax.loglog(t, SigmaDustTotPeak, label="SigmaDust @ peak")
    ax.legend()
    ax.set_xlim(minyear, maxyear)
    ax.set_ylim(1e-3, 1e1)
    filename = "208sds"
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Surface density [g/cm²]")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    parser.add_argument('-t', action="store", dest="title", type=int, default=1, help="Show plot")
    parser.add_argument('-x', action="store", dest="text", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)
