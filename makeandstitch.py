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
    outputDir = localDir + '/figplots/condensedcollages/'

    rows = 4
    columns = 3
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        #dirs = [184, 186, 186, 187, 188, 189, 193, 194, 195, 196, 197, 198]
        # dirs = [184, 185, 186, 187, 188, 189, 190, 191, 192, 196, 197, 198, 199, 200, 201]
        # dirs = [184, 185, 186, 187, 188, 189, 196, 197, 198]
        dirs = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]

    else:
        dirs = [202, 203, 204, 206, 207, 208, 210, 211, 212, 214, 215, 216]
        dirs = [112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123]

    vert = 15
    if args.sdr_stat:
        filename = outputDir + 'sdr_stat_new'
    elif args.mass_stat:
        filename = outputDir + 'mass_stat_final'
    elif args.dist_stat:
        vert = 13
        filename = outputDir + 'dist_stat_new'
    elif args.sdr_mov:
        filename = outputDir + 'sdr_mov_new'
    elif args.mass_mov:
        filename = outputDir + 'mass_mov_final'
    elif args.dist_mov:
        vert = 13
        filename = outputDir + 'dist_mov_new'
    print("Writing...", filename)
    n = columns*rows
    fig, ax = plt.subplots(rows, columns, sharex=True, sharey='row', figsize=(12, vert))
    i = 0
    for r in range(rows):
        for col in range(columns):
            print(dirs[i])
            w.datadir = getDataDir(dirs[i])
            # Surface density
            if args.sdr_stat | args.sdr_mov:
                try:
                    data = w.read.output(150)
                except:
                    try:
                        data = w.read.output(144)
                    except:
                        i = i+1
                        continue
                Nt = np.size(w.read.sequence('t'))
                R = data.grid.r / c.au  # Radial grid cell centers [cm]
                SigmaGas = data.gas.Sigma
                SigmaDust = data.dust.Sigma  # Nr x Nm
                SigmaDustTot = np.sum(SigmaDust, axis=1)  # Nr
                d2g = data.dust.eps
                SigmaPlan = data.planetesimals.Sigma
                PlanMass = data.planetesimals.M / M_earth
                rho_gas = data.gas.rho
                rho_dust = data.dust.rho.sum(-1)
                d2g_mid = rho_dust / rho_gas
                center, width, frac, istart, iend = getRingStats(SigmaPlan, w)
                textstr = getText(PlanMass[-1], center, width, frac)
                ax[r, col].loglog(R, SigmaGas, label="Gas", color=sdr_colors[0])
                ax[r, col].loglog(R, SigmaDustTot, label="Dust", color=sdr_colors[2])
                ax[r, col].loglog(R, SigmaPlan, label="Planetesimals", color=sdr_colors[4])
                ax[r, col].loglog(R, d2g_mid, color=sdr_colors[14], ls='-.', label=r"$\rho_d/\rho_g$", linewidth=1.2)
                if col == (columns - 1) and r == 0:
                    ax[r, col].legend(loc='upper right', frameon=True)
                ax[r, col].set_ylim(1.e-6, 1.e3)
                ax[r, col].set_xlim(10, 170)
                ax[r, col].text(0.03, 0.91, textstr, transform=ax[r, col].transAxes)
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
                DustMass = (np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :]) # Nt x Nr
                DustMassTot = np.sum(DustMass, axis=1)
                # print(DustMassTot[-1] / M_earth)
                DustDiskMassEarth = np.sum(DustMass, axis=1) / M_earth
                PlanDiskMassEarth = w.read.sequence('planetesimals.M') / M_earth
                SigmaPlan = w.read.sequence('planetesimals.Sigma')
                rho_gas = w.read.sequence('gas.rho')
                rho_dust = w.read.sequence('dust.rho').sum(-1)
                d2g_mid = rho_dust / rho_gas

                if args.sdr_stat | args.mass_stat | args.dist_stat:
                    [alpha0, amplitude, position] = getJobParams(dirs[i])
                    iguess = np.argmin(abs(R - position))
                else:
                    iguess = np.argmin(abs(R - 90))

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
                print(position, initialRightDustMass)

                center, width, frac, ist, ien = getRingStats(SigmaPlan[-1], w)
                textstr = getText(PlanDiskMassEarth[-1], center, width, frac, justMass=True, initialExtDust=initialRightDustMass)
                ax[r, col].set_xlim(1.01e4, 1e7)
                # ax[r, col].minorticks_on()
                ax[r, col].set_yticks([1, 10, 100])
                ax[r, col].set_ylim(1e0, 2e5)
                lns1 = ax[r, col].loglog(t, GasDiskMassEarth, label="Gas", color=sdr_colors[2])
                ax2 = ax[r, col].twinx()
                lns5 = ax2.semilogx(t, d2g_mid_at_peak, '-.', linewidth=1.2, color=sdr_colors[14], zorder=2, label=r"$\rho_d/\rho_g(r_{\rm peak})$")
                lns2 = ax[r, col].loglog(t, DustDiskMassEarth, label="Dust", color=sdr_colors[0])
                lns3 = ax[r, col].loglog(t, RingDiskMass, ls='--', label="Ring Dust", color=sdr_colors[1])
                lns4 = ax[r, col].loglog(t, PlanDiskMassEarth, label="Planetesimals", color=sdr_colors[4])
                t = ax2.text(0.04, 0.9, textstr, transform=ax[r, col].transAxes, zorder=1000)
                t.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='white'))
                ax2.set_ylim(0, 1)
                ax2.set_yticks([])
                ax2.set_yticklabels([])

                # Legend
                lns = lns1 + lns2 + lns3 + lns4 + lns5
                labs = [l.get_label() for l in lns]
                if col == (columns - 1) and r == 0:
                    l = ax[r, col].legend(lns, labs, loc='upper right', frameon=True, framealpha=1.)
                    l.set_zorder(2000000)

            elif args.dist_stat | args.dist_mov:
                try:
                    data = w.read.output(150)
                except:
                    try:
                        data = w.read.output(144)
                    except:
                        i = i+1
                        continue
                t = w.read.sequence('t') / c.year
                Nt = t.shape[0]
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
                ax[r, col].text(0.03, 0.91, textstr, color="white", transform=ax[r, col].transAxes, weight='bold', fontweight='bold')
                pltcmap = ax[r, col].contourf(R / c.au, m, np.log10(sigmaD.T), levels=np.linspace(sd_max - 6, sd_max, 7),
                                      cmap="magma", extend="both")
                if i < 6:
                    mmax = 10
                else:
                    mmax = 1e8
                ax[r, col].contour(R / c.au, m, St.T, levels=[1.], colors="white", linewidths=2)
                ax[r, col].contour(R / c.au, m, (St - St_dr[..., None]).T, levels=[0.], colors="C2", linewidths=1.5)
                ax[r, col].contour(R / c.au, m, (St - St_fr[..., None]).T, levels=[0.], colors="C0", linewidths=1.5)
                ax[r, col].set_xlim(13, 170)
                ax[r, col].set_ylim(m[0], mmax)
                ax[r, col].set_xscale("log")
                ax[r, col].set_yscale("log")
                ax[r, col].tick_params(axis='x', labelcolor='black', color='white', which='both')
                ax[r, col].tick_params(axis='y', labelcolor='black', color='white', which='both')
                ax[r, col].xaxis.set_minor_formatter(ScalarFormatter())
                ax[r, col].xaxis.set_major_formatter(ScalarFormatter())
                ax[r, col].xaxis.set_minor_formatter(("%.0f"))
                ax[r, col].set_xticks([30, 60, 90, 150])

                for spine in ax[r, col].spines.values():
                       spine.set_edgecolor('white')
            i += 1

    # Saving figure
    mid = 0.52
    sep = 0.31
    pos_height = 0.97
    title_height = 0.985
    plt.subplots_adjust(left=0.055, bottom=0.03, right=0.98, top=0.95, wspace=0, hspace=0.1)
    if args.dist_stat | args.dist_mov:
        mid = mid - 0.09
        sep = sep - 0.05
    # Initial position or velocity information
    if args.sdr_stat | args.mass_stat | args.dist_stat:
        fig.text(mid - sep, pos_height, r'$r_{\rm g}$ = 30 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$r_{\rm g}$ = 60 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid + sep, pos_height, r'$r_{\rm g}$ = 90 au', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, title_height, "Stationary gap evolved for 10 Myr", ha='center', va='center', fontsize=fontsize)
    else:
        fig.text(mid-sep, pos_height, r'$f = 10\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, pos_height, r'$f = 30\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid+sep, pos_height, r'$f = 100\%$', ha='center', va='center', fontsize=fontsize)
        fig.text(mid, title_height, "Migrating gap initially at 90 au evolved for 10 Myr", ha='center', va='center',
                 fontsize=fontsize)
    fig.text(mid, 0.956, r'$\alpha = 10^{-3}$, A = 3', ha='center', va='center', fontsize=fontsize)
    fig.text(mid, 0.7205, r'$\alpha = 10^{-3}$, A = 10', ha='center', va='center', fontsize=fontsize)
    fig.text(mid, 0.4855, r'$\alpha = 10^{-4}$, A = 3', ha='center', va='center', fontsize=fontsize)
    fig.text(mid, 0.25, r'$\alpha = 10^{-4}$, A = 10', ha='center', va='center', fontsize=fontsize)
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
    plt.savefig(filename+'.png', dpi=300, bbox_inches=p["bbox"], pad_inches=0.02)
    plt.savefig(filename+'.pdf', dpi=300, bbox_inches=p["bbox"], pad_inches=0.02)
    #plt.show()


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

