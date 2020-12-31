import dustpy
from dustpy import plot
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams

import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.gridspec as gridspec
from scipy.signal import argrelextrema
from scipy.interpolate import splrep, sproot, splev
from matplotlib.ticker import ScalarFormatter, StrMethodFormatter, NullFormatter, MultipleLocator, FormatStrFormatter, \
    AutoMinorLocator
import matplotlib.ticker as ticker


M_earth = 5.9722e24 * 1e3  # [g]
localDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'


def main(args):
    z = args.z
    w.datadir = localDir + '/sims/' + str(z)
    outputDir = localDir + '/simplots/'

    # Get basic data from files
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    tMyr = t / 1e6
    tMyrEnd = tMyr[Nt - 1]
    print("tMyrEnd = ", tMyrEnd)
    d2g = w.read.sequence('dust.eps')
    rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]
    m = w.read.sequence('grid.m')  # Mass grid field [g]
    Nm = m.shape[1]  # Number of mass bins
    A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]  # Grid constant
    dm = 2. * (A - 1.) / (A + 1.)  # mass bin width
    r = w.read.sequence('grid.r')  # Radial grid cell centers [cm]
    R = r / c.au  # Radial grid cell centers [AU]
    Nr = R.shape[1]

    SigmaDust = w.read.sequence('dust.Sigma')
    SigmaDustTot = np.sum(SigmaDust, axis=2)
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaPlan = w.read.sequence('planetesimals.Sigma')

    xmin = 10
    xmax = 150
    if z == 87:
        epochs = [0, 2, 4, 5]
        xmax = 100
    elif z == 102:
        epochs = [7, 20, 25, 26]
        xmin = 30

    # plot every epoch
    plotAll = 1
    if plotAll:
        epochs = 0

    ############################# PLOTTING ################################################
    # Choose number of epochs
    cols = 2
    fontsize = 10
    if epochs == 0:
        n = Nt - 1  # To view all make this Nt-1
        gap = 1  # Number of epochs between graphed epochs
        start = 0
        epochs = (np.linspace(start, start + n * gap, n, dtype=int))
        cols = 5
        fontsize = 8
    else:
        n = len(epochs)
        start = epochs[0]

    rows = int(np.ceil((n / cols)))

    fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)
    # fig.subplots_adjust(hspace = .1, wspace=.1)
    fig.text(0.5, 0.01, 'Distance from star [AU]', ha='center')
    fig.text(0.01, 0.5, 'Surface Density [g/cm²]', va='center', rotation='vertical')
    axs = axs.ravel()

    j = 0
    for i in epochs:
        axs[j].loglog(R[i, ...], SigmaDustTot[i, ...], label="Dust")
        axs[j].loglog(R[i, ...], SigmaGas[i, ...], label="Gas")
        axs[j].loglog(R[i, ...], SigmaPlan[i, ...], label="Plan")
        axs[j].set_title("{:.2f} Myrs".format(tMyr[i]), fontsize=fontsize)
        axs[j].set_ylim([1e-4, 1e4])
        axs[j].set_xlim([xmin, xmax])
        axs[j].xaxis.set_major_locator(MultipleLocator(30))
        axs[j].xaxis.set_major_formatter(FormatStrFormatter('%d'))
        axs[j].xaxis.set_minor_locator(MultipleLocator(10))
        axs[j].xaxis.set_minor_formatter(NullFormatter())
        j += 1
    fig.tight_layout()
    plt.savefig(outputDir + 'evol/' + str(z) + '.png', format='png', dpi=300)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=0, help="Simulation number")
    args = parser.parse_args()
    main(args)
