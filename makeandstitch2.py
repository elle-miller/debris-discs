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
from scipy.interpolate import interp1d

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
    outputDir = localDir + '/figplots/condensedcollages/'

    rows = 4
    columns = 3
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        dirs = [184, 185, 186, 187, 188, 189, 193, 194, 195, 196, 197, 198]
        dirs = [184, 185, 186, 187, 188, 189, 193, 194, 195, 196, 197, 198]
        # dirs = [184, 185, 186, 187, 188, 189, 196, 197, 198]
    else:
        dirs = [202, 203, 204, 206, 207, 208, 210, 211, 212, 214, 215, 216]

    if args.sdr_stat:
        filename = outputDir + 'sdr_stat'
    elif args.mass_stat:
        filename = outputDir + 'mass_stat4'
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
    fig, ax = plt.subplots(rows, columns, sharex=True, sharey=True, figsize=(12, 7))
    i = 0
    for r in range(rows):
        for col in range(columns):
            print(dirs[i])
            w.datadir = getDataDir(dirs[i])
            # Surface density
            if args.sdr_stat | args.sdr_mov:
                Nt = np.size(w.read.sequence('t'))
                data = w.read.output(Nt-1)
                R = data.grid.r / c.au  # Radial grid cell centers [cm]
                SigmaGas = data.gas.Sigma
                SigmaDust = data.dust.Sigma  # Nr x Nm
                SigmaDustTot = np.sum(SigmaDust, axis=1)  # Nr
                d2g = data.dust.eps
                SigmaPlan = data.planetesimals.Sigma
                PlanMass = data.planetesimals.M / M_earth

                center, width, frac, istart, iend = getRingStats(SigmaPlan, w)
                textstr = getText(PlanMass[-1], center, width, frac)
                ax[r, col].loglog(R, SigmaGas, label="Gas")
                ax[r, col].loglog(R, SigmaDustTot, label="Dust")
                ax[r, col].loglog(R, SigmaPlan, label="Planetesimals")
                ax[r, col].loglog(R, d2g, ls='--', label="d2g Ratio")
                if col == (columns - 1) and r == 0:
                    ax[r, col].legend(loc='upper right', frameon=True)
                ax[r, col].set_ylim(1.e-6, 1.e3)
                ax[r, col].set_xlim(10, 170)
                ax[r, col].text(0.04, 0.9, textstr, transform=ax[r, col].transAxes)
                ax[r, col].xaxis.set_minor_formatter(ScalarFormatter())
                ax[r, col].xaxis.set_minor_formatter(("%.0f"))
                for axis in [ax[r, col].xaxis]:
                    axis.set_major_formatter(ScalarFormatter())
                    ax[r, col].set_xticks([30, 60, 90, 120])

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
                ax[r, col].loglog(t, GasDiskMassEarth, label="Gas", color="C0")
                ax[r, col].loglog(t, DustDiskMassEarth, label="Dust", color="C1")
                ax[r, col].loglog(t, RingDiskMass, ls='--', label="Ring Dust", color="C4")
                ax[r, col].loglog(t, PlanDiskMassEarth, label="Planetesimals", color="C2")
                if col == (columns - 1) and r == 0:
                    ax[r, col].legend(loc='upper right', frameon=True)
                ax[r, col].set_xlim(1.01e4, 1e7)
                ax[r, col].set_ylim(1e0, 1e3)
                ax[r, col].text(0.02, 0.88, textstr, transform=ax[r, col].transAxes)
                ax[r, col].tick_params(axis='both', which='both')
                #ax2 = ax[r, col].twinx()  # instantiate a second axes that shares the same x-axis
                #color = 'tab:gray'
                #ax2.set_ylabel('Midplane d2g ratio at peak', color=color, rotation=-90)
                #ax2.loglog(t, d2g_mid_at_peak, '-.', linewidth=0.5, color=color)
                #ax2.tick_params(axis='y', labelcolor=color)

            elif args.dist_stat | args.dist_mov:
                t = w.read.sequence('t') / c.year
                Nt = t.shape[0]
                data = w.read.output(Nt - 1)
                PlanMassEarth = data.planetesimals.M[-1] / M_earth
                m = data.grid.m
                R = data.grid.r
                ri = data.grid.ri
                SigmaD = data.dust.Sigma
                SigmaG = data.gas.Sigma
                SigmaPlan = data.planetesimals.Sigma
                eps = data.dust.eps
                cs = data.gas.cs
                delta = data.dust.delta.turb
                OmegaK = data.grid.OmegaK
                St = data.dust.St
                vK = OmegaK * R
                vFrag = data.dust.v.frag
                a = np.mean(m[..., 1:] / m[..., :-1], axis=-1)
                dm = 2. * (a - 1.) / (a + 1.)
                sigmaD = SigmaD[..., :] / dm
                b = vFrag ** 2 / (delta * cs ** 2)
                with np.warnings.catch_warnings():
                    np.warnings.filterwarnings(
                        'ignore',
                        r'invalid value encountered in sqrt')
                    St_fr = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b ** 2))
                p = SigmaG * OmegaK * cs
                _f = interp1d(np.log10(R), np.log10(p), fill_value='extrapolate')
                pi = 10. ** _f(np.log10(ri))
                gamma = np.abs(R / p * np.diff(pi) / np.diff(ri))
                St_dr = eps / gamma * (vK / cs) ** 2
                sd_max = np.ceil(np.log10(sigmaD.max()))
                sg_max = np.ceil(np.log10(SigmaG.max()))
                center, width, frac, istart, iend = getRingStats(SigmaPlan, w)
                textstr = getText(PlanMassEarth, center, width, frac, justMass=True)
                ax[r, col].text(0.04, 0.9, textstr, color="white", transform=ax[r, col].transAxes, fontweight='bold')
                pltcmap = ax[r, col].contourf(R / c.au, m, np.log10(sigmaD.T), levels=np.linspace(sd_max - 6, sd_max, 7),
                                      cmap="magma", extend="both")
                ax[r, col].contour(R / c.au, m, St.T, levels=[1.], colors="white", linewidths=2)
                ax[r, col].contour(R / c.au, m, (St - St_dr[..., None]).T, levels=[0.], colors="C2", linewidths=1.5)
                ax[r, col].contour(R / c.au, m, (St - St_fr[..., None]).T, levels=[0.], colors="C0", linewidths=1.5)
                ax[r, col].axhline(m[-1], color="#AAAAAA", lw=1, ls="--")
                ax[r, col].axvline(R[-1] / c.au, color="#AAAAAA", lw=1, ls="--")
                ax[r, col].set_xlim(R[0] / c.au, 300)
                ax[r, col].set_ylim(m[0], m[-1])
                ax[r, col].set_xscale("log")
                ax[r, col].set_yscale("log")
                ax[r, col].tick_params(axis='x')
                ax[r, col].tick_params(axis='y')
                ax[r, col].xaxis.set_minor_formatter(ScalarFormatter())
                ax[r, col].xaxis.set_major_formatter(ScalarFormatter())
                ax[r, col].xaxis.set_minor_formatter(("%.0f"))
                ax[r, col].set_xticks([30, 60, 90, 150])
            i += 1

    # Saving figure
    mid = 0.50
    far = 0.96
    sep = 0.30
    pos_height = 0.97
    title_height = 0.985
    plt.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.95, wspace=0, hspace=0)
    #plt.subplots_adjust(left=0.05, bottom=0.07, right=0.93, top=0.95, wspace=0, hspace=0) # d2g
    d2g = 0
    if args.dist_stat | args.dist_mov:
        mid = mid - 0.09
        sep = sep - 0.05
    # Initial position or velocity information
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        fig.text(mid - sep, pos_height, r'$r_{\rm p}$ = 30 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$r_{\rm p}$ = 60 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid + sep, pos_height, r'$r_{\rm p}$ = 90 au', ha='center', va='center', fontsize=fontsize)
        # fig.text(mid, title_height, "Stationary pressure trap evolved for 10 Myr", ha='center', va='center', fontsize=fontsize)
    else:
        fig.text(mid-sep, pos_height, r'$f = 10\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$f = 30\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid+sep, pos_height, r'$f = 100\%$', ha='center', va='center', fontsize=fontsize)
        #fig.text(mid, title_height, "Migrating pressure trap initially at 90 au evolved for 10 Myr", ha='center', va='center',
              #   fontsize=fontsize)
    fig.text(far+d2g, 0.84, 'A = 3', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.62, 'A = 10', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.39, 'A = 3', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+d2g, 0.18, 'A = 10', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+0.02+d2g, 0.29, r'$\alpha = 10^{-4}$', ha='center', va='center', rotation=-90, fontsize=fontsize)
    fig.text(far+0.02+d2g, 0.73, r'$\alpha = 10^{-3}$', ha='center', va='center', rotation=-90, fontsize=fontsize)
    if args.sdr_stat | args.sdr_mov:
        fig.text(0, 0.5, 'Surface density [g/cm²]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
    elif args.mass_stat | args.mass_mov:
        fig.text(0, 0.5, r'Mass [$M_\oplus$]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Time [yr]', ha='center', va='center', fontsize=fontsize)
    elif args.dist_stat | args.dist_mov:
        plt.subplots_adjust(left=0.055, bottom=0.03, right=0.85, top=0.95, wspace=0, hspace=0.1)
        fig.text(0, 0.5, 'Particle mass [g]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
        fig.text(mid, 0, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
        cbarcmap = plt.colorbar(pltcmap, ax=ax, fraction=0.05, aspect=70, pad=0.01)
        cbarcmap.ax.set_ylabel("$\log\ \sigma$ [g/cm²]")
    p = getPNG(filename)
    plt.savefig(filename+'.png', dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
    plt.savefig(filename+'.eps', dpi=300, bbox_inches=p["bbox"], pad_inches=0.05)
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

