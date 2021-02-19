from dustpy import plot
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import numpy as np
import argparse
import matplotlib.pyplot as plt
from movieBump import movieBump
from dustpy.plot import panel
from os import path

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

# Read all data in the directory
fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(16, 5))
fontsize = 14
labelsize = 14
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=fontsize)

dir = [184, 185, 186]
i = 0
for z in dir:


    outputDir = localDir + '/simplots/'

    # Time data
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    tMyr = t / 1e6
    tMyrEnd = tMyr[Nt - 1]
    print("tMyrEnd = ", tMyrEnd)

    # Mass data
    m = w.read.sequence('grid.m')  # Mass grid field [g]
    Nm = m.shape[1]  # Number of mass bins
    print(Nm)
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
    SigmaDust = w.read.sequence('dust.Sigma')
    SigmaDustTot = np.sum(SigmaDust, axis=2)
    DustDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / c.M_sun
    DustDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / M_earth
    print("Initial dust disc mass (Earths): ", DustDiskMassEarth[0])
    print("Final dust disc mass (Earths): ", DustDiskMassEarth[-1])
    SigmaDustDist = SigmaDust / dm
    particleSize = w.read.sequence('dust.a')  # Particle size field [cm]

    # Gas information
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
    SigmaGasDist = SigmaGas / dm

    # Create strings for plots
    [alpha, amplitude, position] = getJobParams(z)
    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
    textstr = "Plan Disc Mass: " + str(ptot) + " Earths"
    titlestr = '$r_p$=' + str(position) +' au'

    ax[i].loglog(R[-1, ...], SigmaDustTot[-1, ...], label="Dust")
    ax[i].loglog(R[-1, ...], SigmaGas[-1, ...], label="Gas")
    ax[i].loglog(R[-1, ...], SigmaPlan[-1, ...], label="Planetesimals")
    ax[i].loglog(R[-1, ...], d2g[-1, ...], label="d2g Ratio")
    ax[i].set_ylim(1.e-5, 1.e2)
    ax[i].set_xlim(12, 140)
    #ax[i].set_xlabel("Distance from star [AU]")
    #ax[i].set_ylabel("Surface Density [g/cm²]")
    #ax[i].legend()
    ax[i].set_title(titlestr, font='serif', fontsize=fontsize-2)
    #ax[i].text(0.05, 0.9, textstr, transform=ax[i].transAxes, fontsize=10)
    i = i + 1

fig.text(0.5, 0.04, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
fig.text(0.0, 0.5, 'Surface density [g/cm²]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
fig.tight_layout()
filename = localDir + '/figplots/statcompare'
plt.savefig(filename+'.png', format='png', bbox_inches="tight")
plt.savefig(filename+'.eps', format='eps', bbox_inches="tight")
plt.show()