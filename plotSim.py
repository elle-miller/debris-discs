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


# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

colWidth = 240  # pt
fontsize = 14
labelsize = 14
length = 3
width = 0.5
fig_width_pt = colWidth  # Get this from LaTeX using \showthe\columnwidth


def main(args):
    plt.rcParams.update({'font.size': fontsize})
    plt.style.use('seaborn-paper')
    plt.rc('text', usetex=True)
    plt.rc('xtick', labelsize=labelsize)
    plt.rc('ytick', labelsize=labelsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rcParams["legend.labelspacing"] = 0.3

    # Read all data in the directory
    z = args.z
    print("Sim #", z)
    w.datadir = getDataDir(z)
    outputDir = localDir + '/simplots/'

    if args.plotSDMovie:
        movieBump(z, w.datadir, outputDir + 'movies/sdr/', sd=True)
        return
    if args.plotDistMovie:
        movieBump(z, outputDir + 'movies/sdr/', sd=False)
        return

    # Time data
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    tMyr = t / 1e6
    tMyrEnd = tMyr[Nt - 1]
    print("tMyrEnd = ", tMyrEnd)

    # Mass data
    m = w.read.sequence('grid.m')  # Mass grid field [g]
    Nm = m.shape[1]  # Number of mass bins
    A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]  # Grid constant
    dm = 2. * (A - 1.) / (A + 1.)  # mass bin width

    # Radial data
    r = w.read.sequence('grid.r')  # Radial grid cell centers [cm]
    R = r / c.au  # Radial grid cell centers [AU]
    Nr = R.shape[1]
    rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]

    # Misc data
    alpha = w.read.sequence('gas.alpha')
    d2g = w.read.sequence('dust.eps')

    # Planetesimal information
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
    PlanMass = w.read.sequence('planetesimals.M') / M_earth
    PlanDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / c.M_sun
    PlanDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / M_earth
    print("Mass of final planetesimal disc mass in Earth masses: %.10f" % PlanDiskMassEarth[-1])

    # Dust information
    SigmaDust = w.read.sequence('dust.Sigma') # Nt x Nr x Nm
    SigmaDustTot = np.sum(SigmaDust, axis=2)  # Nt x Nr
    DustMass = (np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :])
    DustDiskMass = np.sum(DustMass, axis=1) / c.M_sun
    DustDiskMassEarth = np.sum(DustMass, axis=1) / M_earth
    particleSize = w.read.sequence('dust.a')  # Particle size field [cm]
    # print("Initial dust disc mass (Earths): ", DustDiskMassEarth[0])
    print("Final dust disc mass (Earths): ", DustDiskMassEarth[-1])

    # Loop through radial grid, setting to 0 on either side of bump
    LeftDustTot = SigmaDustTot[-1].copy()
    RightDustTot = SigmaDustTot[-1].copy()
    iBumpPeakPos = np.argmax(DustMass[-1])
    for i in range(Nr):
        if i <= iBumpPeakPos:
            RightDustTot[i] = 0
        elif i > iBumpPeakPos:
            LeftDustTot[i] = 0
    leftDustDiskMass = np.sum(np.pi * (rInt[-1, 1:] ** 2. - rInt[-1, :-1] ** 2.) * LeftDustTot[:]) / M_earth
    rightDustDiskMass = np.sum(np.pi * (rInt[-1, 1:] ** 2. - rInt[-1, :-1] ** 2.) * RightDustTot[:]) / M_earth
    print("Dust left of bump (Earths): ", leftDustDiskMass)
    print("Dust right of bump (Earths): ", rightDustDiskMass)
    # fig, ax = plt.subplots()
    # ax.plot(R[-1], DustMass[-1] / M_earth)
    # ax.vlines(R[-1, iBumpPeakPos], 0, 0.1, 'r')
    # plt.show()

    # Gas information
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
    SigmaGasDist = SigmaGas / dm

    center, width, frac, istart, iend = getRingStats(SigmaPlan[-1], w)
    RingDiskMass = getRingMass(SigmaDustTot.copy()[-1], istart, iend, w)

    textstr = getText(PlanDiskMassEarth[-1], center, width, frac)
    titlestr = getTitle(z, w)

    # Plot the surface density of dust and gas vs the distance from the star
    if args.plotSDR:
        fig, ax = plt.subplots()
        ax.loglog(R[-1, ...], SigmaGas[-1, ...], label="Gas")
        ax.loglog(R[-1, ...], SigmaDustTot[-1, ...], label="Dust")
        ax.loglog(R[-1, ...], SigmaPlan[-1, ...], label="Plan")
        ax.loglog(R[-1, ...], d2g[-1, ...], ls='--', label="d2g")
        ax.set_ylim(1.e-6, 1.e4)
        ax.set_xlim(12, 200)
        ax.set_xlabel("Distance from star [au]")
        ax.set_ylabel("Surface density [g/cm²]")
        ax.legend(loc='upper right', fontsize=fontsize, frameon=True)
        ax.set_title(titlestr, fontdict={'fontsize': fontsize})
        ax.text(0.04, 0.83, textstr, transform=ax.transAxes, fontsize=fontsize)
        ax.xaxis.set_minor_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(("%.0f"))
        for axis in [ax.xaxis]:
            axis.set_major_formatter(ScalarFormatter())
            ax.set_xticks([30, 60, 90, 120])
        fig.tight_layout()
        filename = outputDir + 'sdr/r' + str(args.z) + '.png'
        plt.savefig(filename, format='png')
        plt.show()

    # Time evolution of gas and dust disk mass
    if args.plotMass:
        fig02, ax02 = plt.subplots()
        ax02.loglog(t, GasDiskMassEarth, label="Gas", color="C1")
        ax02.loglog(t, DustDiskMassEarth, label="Dust", color="C0")
        ax02.loglog(t, RingDiskMass * c.M_sun / M_earth, ls='--', label="Ring Dust", color="C4")
        ax02.loglog(t, PlanDiskMassEarth, label="Planetesimals", color="C3")
        xlim0 = t[min(1, len(t) - 1)]
        xlim1 = t[-1]
        ax02.set_xlim(xlim0, xlim1)
        ax02.set_ylim(1e0, 3e5)
        ax02.legend(loc='lower left')
        ax02.lineTime = ax02.axvline(t[0], color="C7", zorder=-1, lw=1)
        titlestr = r"$\alpha$" + "={a}, A={A}, r$_p$={p}au".format(a=alpha0, A=amplitude, p=position)
        ax02.set_title(titlestr)
        ax02.set_xlabel("Time [yr]")
        ax02.set_ylabel("Mass [M$_\oplus$]")
        ax02.grid(b=False)
        filename = outputDir + 'mass/m' + str(args.z)
        plt.savefig(filename + '.png', format='png', dpi=600)
        plt.savefig(filename + '.eps', format='eps', dpi=600, bbox_inches="tight", pad_inches=0)
        # plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-g', action="store", dest="getRingStats", type=int, default=0, help="Calculate ring stuff")
    parser.add_argument('-a', action="store", dest="plotAll", type=int, default=0, help="Plot all")
    parser.add_argument('-m', action="store", dest="plotMass", type=int, default=0, help="Plot masses over time")
    parser.add_argument('-r', action="store", dest="plotSDR", type=int, default=0, help="Plot sigma over radius")
    parser.add_argument('-s', action="store", dest="plotSDMovie", type=int, default=0, help="Plot a film of colour map")
    parser.add_argument('-d', action="store", dest="plotDistMovie", type=int, default=0, help="Plot plan ring over time")
    arguments = parser.parse_args()
    main(arguments)

# import dustpy
# from dustpy.plotting.plot import readFilesFromDir
# from dustpy.plotting.plot import getSequence
# from dustpy.plotting import plot as specialPlot
# from dustpy.sim.constants import AU
# from dustpy.sim import constants as c
# import numpy as np
# import numpy.ma as ma
# # import numpy.amax as amax
# import matplotlib.pyplot as plt
# from dustpy.plotting.movie import movie
# import argparse
# import matplotlib.gridspec as gridspec
# from scipy.signal import argrelextrema
# from scipy.interpolate import splrep, sproot, splev
# from jobinfo import getJobParams
# from movieBump import movieBump
#
#
# class MultiplePeaks(Exception): pass
#
#
# class NoPeaksFound(Exception): pass
#
#
# M_earth = 5.9722e24 * 1e3  # grams
#
# # Figure stuff
# colWidth = 244  # pt
# fontsize = 12
# labelsize = 14
# length = 3
# width = 0.5
# fig_width_pt = colWidth  # Get this from LaTeX using \showthe\columnwidth
# plt.rcParams.update({'font.size': fontsize})
# plt.rc('axes', labelsize=labelsize)
# plt.rc('xtick', labelsize=labelsize)
# plt.rc('ytick', labelsize=labelsize)
#
#
# def main(args):
#     # Plot flags
#     plotAll = args.plotAll
#     plotMovie = args.plotMovie
#     plotColourS = args.plotColourS
#     plotColourM = 0
#     plotSDR = args.plotSDR
#     plotSDT = args.plotSDT
#     plotMass = args.plotMass
#     plotMoving = args.plotMoving
#     plotPlan = args.plotPlan
#     z = args.outputDirNo
#     outputDir = 'output%s/' % z
#     # if z >300:
#     #     print("hello")
#     # elif z >= 290 | z == 258:
#     #     [alpha, amplitude, velocity] = getJobParams(z)
#     # elif z >= 260:
#     #     [alpha, amplitude, position] = getJobParams(z)
#     # elif z >= 230:
#     #     [alpha, amplitude, velocity] = getJobParams(z)
#     # else:
#     #     [alpha, amplitude, position] = getJobParams(z)
#
#     # Use their plotting function
#     if plotAll:
#         specialPlot(outputDir)
#
#     # if plotMovie:
#     #     fps = 100
#     #     ylim = (1e-12, 1e6)
#     #     showst1 = True
#     #     movie(outputDir, fps=fps, showst1=showst1)
#
#     # Read the files
#     # files = readFilesFromDir(dir="output/",files="*.hdf5")
#     files = readFilesFromDir(dir=outputDir, files="*.hdf5")
#
#     # Extract the data fields from this files, by inputting the field specifier and the output files
#     t = getSequence("t", files) / c.yr
#     tMyr = t / 1e6
#     Nt = t.shape[0]
#     tMyrEnd = tMyr[Nt - 1]
#     d2g = getSequence("dust/dust2gasRatio", files)
#     rInt = getSequence("grid/rInt", files)
#
#     # Read the mass grid m,  has Nt rows and Nm columns
#     #m = getSequence("grid/m", files)
#     #Nm = m.shape[1]
#
#     # Compute the grid constant and the mass bin width
#     #A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]
#     #dm = 2. * (A - 1.) / (A + 1.)
#
#     # Compute the radial grid, convert to AU
#     r = getSequence("grid/r", files)
#     R = r / c.AU
#     Nr = R.shape[1]
#
#     St = getSequence("dust/St", files)
#     deltaTurb = getSequence("dust/deltaTurb", files)
#     if deltaTurb is None:
#         deltaTurb = getSequence("gas/alpha", files)
#     cs = getSequence("gas/cs", files)
#     SigmaGas = getSequence("gas/Sigma", files)
#     omega = getSequence("grid/OmegaK", files)
#     rhos = getSequence("dust/rhoBulk", files).mean(-1)
#     vf = getSequence("dust/vFrag", files)
#     vk = omega * r
#     # calculate the fragmentation limit
#     b = vf ** 2 / (deltaTurb * cs ** 2)
#     with np.warnings.catch_warnings():
#         np.warnings.filterwarnings(
#             'ignore',
#             r'invalid value encountered in sqrt')
#         y_fr0 = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b ** 2))
#
#     # calculate the drift limit
#     p = SigmaGas * omega * cs
#     from scipy.interpolate import interp1d
#     _f = interp1d(np.log10(r[0]), np.log10(p), fill_value='extrapolate')
#     pInt = 10. ** _f(np.log10(rInt[0]))
#     gamma = np.abs(r / p * np.diff(pInt) / np.diff(rInt))
#     y_dr0 = d2g / gamma * (vk / cs) ** 2
#     # Size
#     y_fr_s = 2 * SigmaGas / (np.pi * rhos) * y_fr0
#     y_dr_s = 2 * SigmaGas / (np.pi * rhos) * y_dr0
#     # Mass
#     y_fr_m = 4. * np.pi / 3 * rhos * y_fr0 ** 3
#     y_dr_m = 4. * np.pi / 3 * rhos * y_dr0 ** 3
#
#     # The DustPy units are such that the mass is integrated over a mass bin. That
#     # means numerically we can sum up the mass dimension to get to dust surface density
#     # SigmaDust and SigmaDustTot has Nt rows and Nr cols
#     # Divide by dm to get the logarithmic density distribution
#     SigmaDust = getSequence("dust/Sigma", files)
#     SigmaDustTot = np.sum(SigmaDust, axis=2)
#     DustDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / c.M_sun
#     print(DustDiskMass[0] * c.M_sun / M_earth)
#     #SigmaDustDist = SigmaDust / dm  # 4 x 100
#     particleSize = getSequence("dust/a", files)  # Nt x Nr = 4 x 100
#
#     # Gas information, SigmaGas_m has Nt rows and Nr columns
#     # SigmaGas = getSequence("gas/Sigma", files)
#     SigmaGasTot = np.sum(SigmaGas, axis=-1)
#     GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
#     #SigmaGasDist = SigmaGas / dm
#
#     # Planetesimal information
#     SigmaPlan = getSequence("dust/SigmaPlan", files)
#     SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
#     PlanDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / c.M_sun
#     PlanDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / M_earth
#     print("Plan disk mass in Earth masses: %.1f" % PlanDiskMassEarth[Nt - 1])
#
#     # Find first instance of planetesimals
#     # for j in range(Nt):
#     #    for i in range(Nr):
#     #        if SigmaPlan[j,i]>1e-99:
#     #            print("bingo")
#     #            print("time = %f" % (t[i]/1e6))
#     #            print("pos = %f" % R[j,i])
#     #            print("d2g = %f" % d2g[j,i])
#     # Initialize
#     ringWidth = np.zeros(Nt)
#     startRing = np.zeros(Nt)
#     endRing = np.zeros(Nt)
#     centerRing = np.zeros(Nt)
#     fracWidth = np.zeros(Nt)
#     ringWidth2 = np.zeros(Nt)
#     startRing2 = np.zeros(Nt)
#     endRing2 = np.zeros(Nt)
#     centerRing2 = np.zeros(Nt)
#     fracWidth2 = np.zeros(Nt)
#
#     # Duss in Ring
#     RingDustTot = SigmaDustTot.copy()
#     Ring1DustTot = SigmaDustTot.copy()
#     Ring2DustTot = SigmaDustTot.copy()
#
#     # For every epoch
#     twoRingsFlag = False  # Bad coding practice, basically once two rings are flagged this forces always two rings
#     formed = False
#     for j in range(Nt):
#
#         # Reset index positions
#         istartRingP = 0
#         iendRingP = 0
#         iStartRing2P = 0
#         iEndRing2P = 0
#         index = 0
#         floorVal = 1.0e-90
#         minVal = min(SigmaPlan[j])  # 1e-100
#
#
#         # Loop through each radial bin, locating index positions of start and end ring
#         numRings = 0
#         beginRing2 = False
#         for i in SigmaPlan[j]:
#             if (i > minVal) & (istartRingP == 0):
#                 if formed is False:
#                     formationTimeIndex = j
#                     formed = True
#                 istartRingP = index
#                 numRings = 1
#             elif (beginRing2 == False) & (numRings == 1) & (i < floorVal) & (istartRingP != 0):
#                 iendRingP = index
#                 beginRing2 = True
#             elif beginRing2 & (i > minVal) & (iStartRing2P == 0):
#                 iStartRing2P = index
#                 numRings = 2
#             elif beginRing2 & (i < floorVal) & (iStartRing2P != 0):
#                 iEndRing2P = index
#                 break
#             index += 1
#
#         startRing[j] = 1e-10
#         endRing[j] = 1e-10
#         startRing2[j] = 1e-10
#         endRing2[j] = 1e-10
#
#         # Convert these indices to actual values
#         if istartRingP != 0:
#             startRing[j] = rInt[j, istartRingP] / c.AU
#             endRing[j] = rInt[j, iendRingP] / c.AU
#
#         if iStartRing2P != 0:
#             startRing2[j] = rInt[j, iStartRing2P] / c.AU
#             endRing2[j] = rInt[j, iEndRing2P] / c.AU
#
#         centerRing[j] = (endRing[j] + startRing[j]) / 2
#         ringWidth[j] = endRing[j] - startRing[j]
#         fracWidth[j] = ringWidth[j] / centerRing[j]
#         centerRing2[j] = (endRing2[j] + startRing2[j]) / 2
#         ringWidth2[j] = endRing2[j] - startRing2[j]
#         fracWidth2[j] = ringWidth2[j] / centerRing2[j]
#
#         # Now set values not in dust ring to 0
#         # fw=fwhm(R[j],SigmaPlan[j])
#         # iStartRingP = (np.abs(R - fw[1])).argmin()
#         # iEndRingP =  (np.abs(R - fw[0])).argmin()
#
#         # These are the index values of dust calculated using different way, FWHM
#         # fw=fwhm(R[j],SigmaDustTot[j])
#         # iStartRingD = (np.abs(R - fw[1])).argmin()
#         # iEndRingD =  (np.abs(R - fw[0])).argmin()
#
#         # For this epoch, loop through all radial bins, and turn OFF dust outside ring indices
#         # print(numRings)
#         # print(istartRingP)
#         # print(iendRingP)
#         # print(iStartRing2P)
#         # print(iEndRing2P)
#         for k in range(Nr):
#             # If there are no rings, turn off everything
#             if (numRings == 0):
#                 RingDustTot[j, :] = 0
#                 Ring1DustTot[j, :] = 0
#                 Ring2DustTot[j, :] = 0
#                 break
#             # If there is only one ring and we aint in it
#             elif (numRings == 1) & ((k not in range(istartRingP, iendRingP + 1))):
#                 RingDustTot[j, k] = 0
#                 Ring1DustTot[j, k] = 0
#                 Ring2DustTot[j, :] = 0
#             # If there are two rings, and we are in first, middle or end section then set to zero
#             elif (numRings == 2) & ((k < istartRingP) | (k in range(iendRingP + 1, iStartRing2P)) | (k > iEndRing2P)):
#                 RingDustTot[j, k] = 0
#                 Ring1DustTot[j, k] = 0
#                 Ring2DustTot[j, k] = 0
#             # If in first ring, set ring2 to zero
#             elif (numRings == 2) & (k in range(istartRingP, iendRingP + 1)):
#                 Ring2DustTot[j, k] = 0
#             elif (numRings == 2) & (k in range(iStartRing2P, iEndRing2P + 1)):
#                 Ring1DustTot[j, k] = 0
#
#     RingDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * RingDustTot[:, :], axis=1) / c.M_sun
#     Ring1DiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * Ring1DustTot[:, :], axis=1) / c.M_sun
#     Ring2DiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * Ring2DustTot[:, :], axis=1) / c.M_sun
#    # print("Time of first formation: %.2f" % tMyr[formationTimeIndex])
#     print("Ring 1 Start: %.1f" % startRing[-1])
#     print("Ring 1 End: %.1f" % endRing[-1])
#     print("Ring 1 Center: %.1f AU" % centerRing[Nt - 1])
#     print("Ring 1 Width: %.1f AU" % ringWidth[Nt - 1])
#     print("Ring 1 Fractional Width: %.2f" % fracWidth[Nt - 1])
#     print("******")
#     for i in range(Nt-1):
#         print("Time = %.2f" % tMyr[i])
#         print("Plan Mass = %.2f" % PlanDiskMassEarth[i])
#
#     if ringWidth2[Nt - 1] > 0:
#         print("Ring 2 Center: %.1f AU" % centerRing2[Nt - 1])
#         print("Ring 2 Width: %.1f AU" % ringWidth2[Nt - 1])
#         print("Ring 2 Fractional Width: %.2f" % fracWidth2[Nt - 1])
#
#     c1 = f"{centerRing[Nt - 1]:.1f}"
#     w1 = f"{ringWidth[Nt - 1]:.1f}"
#     f1 = f"{fracWidth[Nt - 1]:.2f}"
#     c2 = f"{centerRing2[Nt - 1]:.0f}"
#     w2 = f"{ringWidth2[Nt - 1]:.0f}"
#     f2 = f"{fracWidth2[Nt - 1]:.2f}"
#     ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
#
#     if ringWidth[Nt - 1] == 0:
#         textstr = ""
#     elif ringWidth2[Nt - 1] == 0:
#         textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1
#     else:
#         textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1 + "\n" + "2: c=" + c2 + " w=" + w2 + " AU, f=" + f2
#
#     titlestr = ""
#     if ((z > 200) & (z < 260)) | (z >= 290) | (z == 258):
#         plottingVel = True
#         textstr = textstr + "\n" + "Plan Disc Mass: " + str(ptot) + " Earths"
#         titlestr = str(z) + "at " + str(tMyrEnd) +"Myrs" #+ ": " + r"$\alpha$" + "={a}, A={A}, v={v}% @ {t:.2f} Myr".format(a=alpha, A=amplitude,
#                                                                                      #       v=velocity,
#                                                                                      #       t=tMyrEnd)
#     else:
#         plottingVel = False
#         textstr = textstr + "\n" + "Plan Disc Mass: " + str(ptot) + " Earths"
#         titlestr = str(z) #+ ": " + r"$\alpha$" + "={a}, A={A}, $r_p$={p}AU @ {t:.2f} Myr".format(a=alpha, A=amplitude,
#                             #                                                                     p=position,
#                               #                                                                   t=tMyrEnd)
#
#     ############################# PLOTTING ################################################
#

#
#     # Plot the surface density of dust and gas vs the distance from the star
#     if plotSDR:
#         fig, ax = plt.subplots()
#         ax.loglog(R[-1, ...], SigmaDustTot[-1, ...], label="Dust")
#         ax.loglog(R[-1, ...], SigmaGas[-1, ...], label="Gas")
#         ax.loglog(R[-1, ...], SigmaPlan[-1, ...], label="Planetesimals")
#         ax.loglog(R[-1, ...], d2g[-1, ...], label="d2g Ratio")
#         ylim0 = 10. ** np.floor(np.log10(np.min(np.minimum(SigmaDustTot, SigmaGas))))
#         ax.set_ylim(1.e-6, 1.e4)
#         ax.set_xlabel("Distance from star [AU]")
#         ax.set_ylabel("Surface Density [g/cm²]")
#         ax.legend()
#         ax.set_title(titlestr)
#         ax.text(0.05, 0.9, textstr, transform=ax.transAxes, fontsize=10)
#         fig.tight_layout()
#         # plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/r%d.eps' % z, format='eps',
#         # dpi=1200)
#         plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/r%d.png' % z, format='png',
#                     dpi=600)
#         plt.show()
#
#     if plotPlan:
#         fig = plt.figure()
#         gs = gridspec.GridSpec(2, 2)
#         ax2 = fig.add_subplot(gs[0, 0])
#         ax3 = fig.add_subplot(gs[0, 1])
#         ax1 = fig.add_subplot(gs[1, :])
#
#         # ax1.loglog(R[:,...],SigmaPlan[:,...])
#         ax1.set_xlabel("Distance from star [AU]")
#         ax1.set_ylabel("Plan SD [g/cm²]")
#         ax1.set_ylim(1.e-4, 1e0)
#         # ax1.set_xlim(10,6e1)
#         for it in range(Nt - 6, Nt):
#             ax1.loglog(R[0, ...], SigmaPlan[it, ...], label='%d yrs' % t[it])
#             ax1.legend(prop={'size': 8})
#
#         # plot,ax = plt.subplots(1,2)
#         ax2.set_xlabel("Time [Myr]")
#         ax2.set_ylabel("Ring Width [AU]")
#         ax2.plot(t / 1.0e6, ringWidth)
#
#         # p#lot,ax = plt.subplots(1,3)
#         ax3.plot(t / 1.0e6, fracWidth)
#         ax3.plot(t / 1.0e6, np.ones(Nt) * 0.5)
#         ax3.set_ylabel('Fractional Width')
#         ax3.set_xlabel("Time [Myr]")
#         plt.show()
#
#     # Plot the dust density distribution vs mass with different times
#     if plotSDT:
#         plot01, ax01 = plt.subplots()
#         for it in range(Nt):
#             plot01 = ax01.loglog(m[0, ...], SigmaDustDist[it, 0, ...], label='%d yrs' % t[it])
#         ax01.set_xlabel("Mass [g]")
#         ax01.set_ylabel("$\sigma_\mathrm{dust}$ [g/cm$^2$]")
#         ax01.set_xlim(min(m[0]))
#         ax01.legend()
#         plt.show()
#
#     # Calculate the limits of the logarithmic distribution to manage the colorbar
#     if plotColourS:
#         SigmaDustDistMax = np.ceil(np.max(np.log10(SigmaDustDist)))
#         levels = np.linspace(SigmaDustDistMax - 9., SigmaDustDistMax, 10)
#         fig00, ax00 = plt.subplots()
#         plot = ax00.contourf(R[-1][:, None] * np.ones_like(particleSize[-1, ...]), particleSize[-1, ...],
#                              np.log10(SigmaDustDist[-1, ...]), levels=levels, extend="both")
#         plot00_St1 = ax00.contour(R[-1][:, None] * np.ones_like(particleSize[-1, ...]), particleSize[-1, ...],
#                                   St[-1, ...], levels=[1.], colors="#FFFFFF")
#         ax00St1Collections = plot00_St1.collections[:]
#         plot00_Stdr = ax00.contour(R[-1][:, None] * np.ones_like(particleSize[-1, ...]), particleSize[-1, ...],
#                                    (St[-1, ...] - y_dr_s[-1, ..., None]), levels=[0.], colors="C9", linewidths=1.)
#         ax00StdrCollections = plot00_Stdr.collections[:]
#         plot00_Stf = ax00.contour(R[-1][:, None] * np.ones_like(particleSize[-1, ...]), particleSize[-1, ...],
#                                   (St[-1, ...] - y_fr_s[-1, ..., None]), levels=[0.], colors="C6", linewidths=1.)
#         ax00StfCollections = plot00_Stf.collections[:]
#         cbar = plt.colorbar(plot, ax=ax00)
#         ax00.set_ylim(1e-4, 1e2)
#         ax00.set_xscale("log")
#         ax00.set_yscale("log")
#         ax00.set_xlabel("Distance from star [AU]")
#         ax00.set_ylabel("Particle radius [cm]")
#         cbar.ax.set_ylabel("log $\sigma_\mathrm{dust}$ [g/cm²]")
#         plt.show()
#
#     # Plot mass of planetesimals
#     if plotColourM:
#         SigmaDustDistMax = np.ceil(np.max(np.log10(SigmaDustDist)))
#         levels = np.linspace(SigmaDustDistMax - 9., SigmaDustDistMax, 10)
#         fig00, ax00 = plt.subplots()
#         plot = ax00.contourf(R[-1, ...], m[-1, ...], np.log10(SigmaDustDist[-1, ...].T), levels=levels, extend="both")
#         plot00_St1 = ax00.contour(R[-1, ...], m[-1, ...], St[2, ...].T, levels=[1.], colors="#FFFFFF")
#         ax00St1Collections = plot00_St1.collections[:]
#         plot00_Stdr = ax00.contour(R[-1, ...], m[-1, ...], (St[2, ...] - y_dr_s[2, ..., None]).T, levels=[0.],
#                                    colors="C9", linewidths=1.)
#         ax00StdrCollections = plot00_Stdr.collections[:]
#         plot00_Stf = ax00.contour(R[-1, ...], m[-1, ...], (St[2, ...] - y_fr_s[2, ..., None]).T, levels=[0.],
#                                   colors="C6", linewidths=1.)
#         ax00StfCollections = plot00_Stf.collections[:]
#         cbar = plt.colorbar(plot, ax=ax00)
#         ax00.set_xscale("log")
#         ax00.set_yscale("log")
#         ax00.set_xlabel("Distance from star [AU]")
#         ax00.set_ylabel("Mass [g]")
#         cbar.ax.set_ylabel("log $\sigma_\mathrm{dust}$ [g/cm²]")
#         # plt.show()
#
#     # Time evolution of gas and dust disk mass
#     if plotMass:
#
#         fig02, ax02 = plt.subplots()
#         ax02.loglog(t, GasDiskMass * c.M_sun / M_earth, label="Gas", color="C0")
#         ax02.loglog(t, DustDiskMass * c.M_sun / M_earth, label="Dust", color="C4")
#         if numRings == 2:
#             ax02.loglog(t, RingDiskMass * c.M_sun / M_earth, ls='--', label="Total Ring Dust", color="C1")
#             ax02.loglog(t, Ring1DiskMass * c.M_sun / M_earth, ls='-.', label="Ring 1 Dust", color="C3")
#             ax02.loglog(t, Ring2DiskMass * c.M_sun / M_earth, ls=':', label="Ring 2 Dust", color="C5")
#         else:
#             ax02.loglog(t, RingDiskMass * c.M_sun / M_earth, ls='--', label="Ring Dust", color="C1")
#         ax02.loglog(t, PlanDiskMass * c.M_sun / M_earth, label="Planetesimals", color="C2")
#         xlim0 = t[min(1, len(t) - 1)]
#         xlim1 = t[-1]
#         ax02.set_xlim(xlim0, xlim1)
#         # ylim0 = 10. ** np.floor(np.log10(np.min(np.minimum(GasDiskMass, DustDiskMass, PlanDiskMass))))
#         # ylim1 = 10. ** np.ceil(np.log10(np.max(np.maximum(GasDiskMass, DustDiskMass, PlanDiskMass))))
#         ax02.set_ylim(1e0, 3e5)
#         ax02.legend(loc='lower left')
#         ax02.lineTime = ax02.axvline(t[0], color="C7", zorder=-1, lw=1)
#         ax02.set_title(titlestr)
#         ax02.set_xlabel("Time [yr]")
#         ax02.set_ylabel("Mass [M$_\oplus$]")
#         ax02.grid(b=False)
#         plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/earthmass' + str(z) + '.png',
#                     format='png', dpi=600)
#         plt.show()
#
#     if plotMoving:
#         fig, ax = plt.subplots()
#         ax.set_ylim(1.e-3, 1.e2)
#         nsnap = Nt - 1
#         cmap = plt.get_cmap('viridis')
#         # ax.set_xlim(min(r),max(r))
#         for i in range(nsnap):
#             ax.plot(r[i, :] / c.AU, SigmaGas[i, :], color=cmap(i / nsnap), label='%d yrs' % t[i])
#             # ax.plot(r[i,:]/c.AU, SigmaDustTot[i,:], color=cmap(i/nsnap),label= '%d yrs' % t[i])
#             ax.set_xlim(30, 400)
#             ax.set_xlabel("Distance from star [au]")
#             ax.set_ylabel("Gas Surface Density [g/cm²]")
#             ax.set_xscale('log')
#             ax.set_yscale('log')
#             ax.legend()
#         plt.show()
#
#
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-z', action="store", dest="outputDirNo", type=int, default=0, help="Simulation number")
#     parser.add_argument('-a', action="store", dest="plotAll", type=int, default=0, help="Plot all")
#     parser.add_argument('-m', action="store", dest="plotMass", type=int, default=0, help="Plot masses over time")
#     parser.add_argument('-r', action="store", dest="plotSDR", type=int, default=0, help="Plot sigma over radius")
#     parser.add_argument('-t', action="store", dest="plotSDT", type=int, default=0, help="Plot sigma over time")
#     parser.add_argument('-c', action="store", dest="plotColourS", type=int, default=0,
#                         help="Plot colour map against particle size")
#     parser.add_argument('-f', action="store", dest="plotMovie", type=int, default=0, help="Plot a film of colour map")
#     parser.add_argument('-g', action="store", dest="plotMoving", type=int, default=0, help="Plot sigma gas over time")
#     parser.add_argument('-p', action="store", dest="plotPlan", type=int, default=0,
#                         help="Plot planetesimal ring over time")
#     args = parser.parse_args()
#     main(args)
#
#
# # Gas distribution
# # SigmaGasDistMax = np.ceil(np.max(np.log10(SigmaGasDist)))
# # levels = np.linspace(SigmaGasDistMax-9., SigmaGasDistMax, 10)
# # fig3, ax3 = plt.subplots()
# # plot = ax3.contourf(R[-1, ...], m[-1, ...], np.log10(SigmaGasDist[-1, ...].T), levels=levels, extend="both")
# # cbar = plt.colorbar(plot, ax=ax3)
# # ax3.set_xscale("log")
# # ax3.set_yscale("log")
# # ax3.set_xlabel("Distance from star [AU]")
# # ax3.set_ylabel("Mass [g]")
# # cbar.ax.set_ylabel("log $\sigma_\mathrm{gas}$ [g/cm²]")
# # plt.show()
#
# # Time evolution
# # fig, ax = plt.subplots()
# # ax.loglog(t[-1, ...], SigmaDust[-1, ...],label="Dust")
# # ax.loglog(t[-1, ...], SigmaGas_m[-1, ...],label="Gas")
# # ax.set_ylim(1.e-4, 1.e4)
# # x.set_xlabel("Distance from star [AU]")
# # ax.set_ylabel("Surface Density [g/cm²]")
# # ax.legend()
# # plt.show()
#
#
