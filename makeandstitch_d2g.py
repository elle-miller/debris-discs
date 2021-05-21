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
from matplotlib.ticker import ScalarFormatter, LogLocator, NullFormatter
from plottingFunctions import *
from dustpy import hdf5writer as w
from scipy.interpolate import interp1d
from matplotlib.cm import get_cmap

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

sdr_colors = [i for i in get_cmap('tab20').colors]

def main(args):

    # Overlay seaborn's styling with personal adjustments
    plt.style.use('seaborn-paper')
   # plt.style.use('tex')
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

    rows = 4
    columns = 3
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        dirs = [184, 184, 186, 187, 188, 189, 193, 194, 195, 196, 197, 198]
        dirs = [184, 184, 186, 187, 188, 189, 193, 194, 195, 196, 197, 198]
        # dirs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        # dirs = [184, 185, 186, 187, 188, 189, 196, 197, 198]
    else:
        dirs = [202, 203, 204, 206, 207, 208, 210, 211, 212, 214, 215, 216]

    dirs = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
    if args.sdr_stat:
        filename = outputDir + 'sdr_stat'
    elif args.mass_stat:
        filename = outputDir + 'mass_stat_new'
    elif args.dist_stat:
        filename = outputDir + 'dist_stat'
    elif args.sdr_mov:
        filename = outputDir + 'sdr_mov'
    elif args.mass_mov:
        filename = outputDir + 'mass_mov'
    elif args.dist_mov:
        filename = outputDir + 'dist_mov'
    print("Writing...", filename)
    n = columns*rows
    fig, ax = plt.subplots(rows, columns, sharex=True, sharey=True, figsize=(12, 8))
    i = 0
    for r in range(rows):
        for col in range(columns):
            print(dirs[i])
            w.datadir = getDataDir(dirs[i])
            # Surface density
            if args.sdr_stat | args.sdr_mov:
                Nt = np.size(w.read.sequence('t'))

            # Mass plotting
            elif args.mass_stat | args.mass_mov:
                t = w.read.sequence('t') / c.year
                Nt = t.shape[0]
                R = w.read.sequence('grid.r') / c.au  # Radial grid cell centers [cm]
                Nr = R.shape[1]
                rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]
                SigmaGas = w.read.sequence('gas.Sigma')
                GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :],
                                           axis=1) / M_earth
                SigmaDust = w.read.sequence('dust.Sigma')  # Nt x Nr x Nm
                SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
                DustMass = (np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :])
                DustDiskMassEarth = np.sum(DustMass, axis=1) / M_earth
                PlanDiskMassEarth = w.read.sequence('planetesimals.M') / M_earth
                SigmaPlan = w.read.sequence('planetesimals.Sigma')
                rho_gas = w.read.sequence('gas.rho')
                rho_dust = w.read.sequence('dust.rho').sum(-1)
                d2g_mid = rho_dust / rho_gas

                if args.sdr_stat | args.mass_stat | args.dist_stat:
                    [alpha0, amplitude, position] = getJobParams(dirs[i])
                    iguess = np.argmin(abs(R[-1] - position))
                else:
                    iguess = np.argmin(abs(R[-1] - 90))

                # Obtain initial external dust
                it = 0
                initialRightDustTot = SigmaDustTot[it].copy()
                for ir in range(Nr):
                    if ir <= iguess:
                        initialRightDustTot[ir] = 0
                initialRightDustMass = np.sum(
                    np.pi * (rInt[it, 1:] ** 2. - rInt[it, :-1] ** 2.) * initialRightDustTot[:]) / M_earth
                d2g_mid_at_peak = np.zeros(Nt)
                r_peak = np.zeros(Nt)

                RingDiskMass = np.zeros(Nt)
                RingDustTot = SigmaDustTot.copy()
                for k in range(Nt):
                    igap = np.argmin(SigmaGas[k, 0:iguess + 10])
                    iguess = igap
                    ipeak = np.argmax(SigmaDustTot[k, igap:igap + 35]) + igap
                    dist = ipeak - igap
                    istart = int(ipeak - 0.5 * dist)
                    iend = int(ipeak + 0.5 * dist)
                    r_peak[k] = R[k, ipeak]
                    d2g_mid_at_peak[k] = d2g_mid[k, ipeak]
                    if ipeak != igap:
                        RingDiskMass[k] = getRingMass(RingDustTot[k], istart, iend, w)
                center, width, frac, ist, ien = getRingStats(SigmaPlan[-1], w)
                textstr = getText(PlanDiskMassEarth[-1], center, width, frac, justMass=True, initialExtDust=initialRightDustMass)
                ax[r, col].set_xlim(1.01e4, 1e7)
                ax[r, col].set_ylim(1e0, 5e2)
                ax[r, col].minorticks_on()
                ax[r, col].set_yticks([1, 10, 100])

                lns1 = ax[r, col].loglog(t, GasDiskMassEarth, label="Gas", color=sdr_colors[2])
                ax2 = ax[r, col].twinx()
                lns5 = ax2.semilogx(t, d2g_mid_at_peak, '-.', linewidth=0.9, color=sdr_colors[14], zorder=2, label=r"$\rho_d/\rho_g(r_{\rm peak})$")
                lns2 = ax[r, col].loglog(t, DustDiskMassEarth, label="Dust", color=sdr_colors[0])
                lns3 = ax[r, col].loglog(t, RingDiskMass, ls='--', label="Ring Dust", color=sdr_colors[1])
                lns4 = ax[r, col].loglog(t, PlanDiskMassEarth, label="Planetesimals", color=sdr_colors[4])
                t = ax2.text(0.04, 0.07, textstr, transform=ax[r, col].transAxes, zorder=1000)
                t.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='white', pad=0.0))
                ax2.set_ylim(0, 1)
                ax2.set_yticks([])
                ax2.set_yticklabels([])

                # Legend
                lns = lns2 + lns3 + lns4 + lns5
                labs = [l.get_label() for l in lns]
                if col == (columns - 2) and r == 0:
                    l = ax2.legend(lns, labs, loc='upper right', frameon=True, framealpha=0.8)
                    ax2.legend(lns, labs, loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=4, fancybox=False)
                    l.set_zorder(2000000)

                # locmaj = LogLocator(base=10, numticks=3)
                # ax[r, col].yaxis.set_major_locator(locmaj)
                # locmin = LogLocator(base=10.0, subs=(0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10))
                # ax[r, col].yaxis.set_minor_locator(locmin)
                # ax[r, col].yaxis.set_minor_formatter(NullFormatter())
            i += 1

    # Saving figure
    mid = 0.50
    far = 0.96
    sep = 0.30
    pos_height = 0.92
    title_height = 0.985
    plt.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.90, wspace=0, hspace=0)
    #plt.subplots_adjust(left=0.05, bottom=0.07, right=0.93, top=0.95, wspace=0, hspace=0) # d2g
    d2g = 0
    if args.dist_stat | args.dist_mov:
        mid = mid - 0.09
        sep = sep - 0.05
    # Initial position or velocity information
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        fig.text(mid - sep, pos_height, r'$r_{\rm g}$ = 30 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$r_{\rm g}$ = 60 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid + sep, pos_height, r'$r_{\rm g}$ = 90 au', ha='center', va='center', fontsize=fontsize)
        # fig.text(mid, title_height, "Stationary pressure trap evolved for 10 Myr", ha='center', va='center', fontsize=fontsize)
    else:
        fig.text(mid-sep, pos_height, r'$f = 10\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$f = 30\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid+sep, pos_height, r'$f = 100\%$', ha='center', va='center', fontsize=fontsize)
        #fig.text(mid, title_height, "Migrating pressure trap initially at 90 au evolved for 10 Myr", ha='center', va='center',
              #   fontsize=fontsize)
    fig.text(far+d2g, 0.8, 'A = 3', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.58, 'A = 10', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.38, 'A = 3', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.17, 'A = 10', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+0.02+d2g, 0.29, r'$\alpha = 10^{-4}$', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+0.02+d2g, 0.73, r'$\alpha = 10^{-3}$', ha='center', va='center', rotation=-90, fontsize=fontsize)
    if args.sdr_stat | args.sdr_mov:
        fig.text(0, 0.5, 'Surface density [g/cm²]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
    elif args.mass_stat | args.mass_mov:
        fig.text(0, 0.5, r'Mass [M$_\oplus$]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Time [yr]', ha='center', va='center', fontsize=fontsize)
    elif args.dist_stat | args.dist_mov:
        plt.subplots_adjust(left=0.055, bottom=0.03, right=0.85, top=0.95, wspace=0, hspace=0.1)
        fig.text(0, 0.5, 'Particle mass [g]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
        cbarcmap = plt.colorbar(pltcmap, ax=ax, fraction=0.05, aspect=70, pad=0.01)
        cbarcmap.ax.set_ylabel("$\log\ \sigma$ [g/cm²]")
    plt.savefig(filename + '.pdf', dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.savefig(filename+'.png', dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.savefig(filename+'.eps', dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-1', action="store", dest="sdr_stat", type=int, default=0, help="Show plot")
    parser.add_argument('-2', action="store", dest="mass_stat", type=int, default=0, help="Show plot")
    parser.add_argument('-3', action="store", dest="dist_stat", type=int, default=0, help="Show plot")
    parser.add_argument('-4', action="store", dest="sdr_mov", type=int, default=0, help="Show plot")
    parser.add_argument('-5', action="store", dest="mass_mov", type=int, default=0, help="Show plot")
    parser.add_argument('-6', action="store", dest="dist_mov", type=int, default=0, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)

