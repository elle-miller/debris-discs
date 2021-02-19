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
fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(5, 12))
fontsize = 14
labelsize = 14
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=fontsize)
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=fontsize)
firstendring = 115
dir = [206, 207, 208]
i = 0
for z in dir:

    if path.exists(localDir + '/sims/' + str(z)):
        dataDir = localDir + '/sims/' + str(z)
    elif path.exists(localDirNew + '/' + str(z)):
        dataDir = localDirNew + '/' + str(z)
    else:
        print("Output directory not found")
        exit()
    w.datadir = dataDir
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
    # print("Initial dust disc mass (Earths): ", DustDiskMassEarth[0])
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
    titlestr = r'$\mathit{f} = $' + str(position) + '%'

    # Duss in Ring
    RingDustTot = SigmaDustTot.copy()
    j = Nt-1

    # Reset index positions
    istartRingP = 0
    iendRingP = 0
    index = 0
    minVal = 1e-3

    # Loop through each radial bin, locating index positions of start and end ring
    for k in SigmaPlan[j]:
        if (k > minVal) & (istartRingP == 0):
            istartRingP = index
        elif (k < minVal) & (istartRingP != 0):
            iendRingP = index
            break
        index += 1

    startRing = 1e-10
    endRing = 1e-10

    # Convert these indices to actual values
    if istartRingP != 0:
        startRing = rInt[j, istartRingP] / c.au
        endRing = rInt[j, iendRingP] / c.au

    centerRing = (endRing + startRing) / 2
    ringWidth = endRing - startRing
    fracWidth = ringWidth / centerRing

    c1 = f"{centerRing:.0f}"
    w1 = f"{ringWidth:.0f}"
    f1 = f"{fracWidth:.2f}"
    if ringWidth != 0:
        textstr1 = "\n" + "c: " + c1 + "au, w: " + w1 + "au, f: " + f1
    else:
        textstr1 = ""

    # Create strings for plots
    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
    textstr = "Planetesimal disc: " + str(ptot) + r" M$_{\oplus}$" + textstr1

    ax[i].loglog(R[-1], SigmaDustTot[-1], label="Dust")
    ax[i].loglog(R[-1], SigmaGas[-1], label="Gas")
    ax[i].loglog(R[-1, ...], SigmaPlan[-1, ...], label="Planetesimals")
    ax[i].loglog(R[-1, ...], d2g[-1, ...], label="d2g Ratio")
    ax[i].vlines(firstendring, 1e-5, 1e3, 'r', '--')
    ax[i].set_ylim(1.e-5, 1.e3)
    ax[i].set_xlim(12, 140)
    ax[i].set_title(titlestr, font='serif', fontsize=fontsize-2)
    ax[i].text(0.05, 0.83, textstr, transform=ax[i].transAxes, fontsize=12)
    i = i + 1

fig.text(0.5, 0.01, 'Distance from star [au]', ha='center', va='center', fontsize=fontsize)
fig.text(0.02, 0.5, 'Surface density [g/cm²]', ha='center', va='center', rotation='vertical', fontsize=fontsize)
fig.tight_layout()
filename = localDir + '/figplots/movingcompare3'
plt.savefig(filename+'.png', format='png', bbox_inches="tight", pad_inches=0)
plt.savefig(filename+'.eps', format='eps', bbox_inches="tight", pad_inches=0)
plt.show()