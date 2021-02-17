import dustpy.constants as c
from dustpy.simulation import Simulation

import numpy as np
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from simframe.io.writers import hdf5writer as writer
from os import path, getcwd
from plottingFunctions import *
from matplotlib.ticker import ScalarFormatter

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

width_inches = 6
golden_ratio = 3 / 4
height_inches = golden_ratio * width_inches
fontsize = 14


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

    t = w.read.sequence('t') / c.year
    tMyrEnd = t[-1] * 1e-6
    Nt = t.shape[0]

    data = w.read.output(Nt - 1)
    filename = "data"
    extension = "hdf5"
    PlanMassEarth = data.planetesimals.M[-1] / M_earth

    Nm = data.grid.Nm
    Nr = data.grid.Nr
    m = data.grid.m
    r = data.grid.r
    ri = data.grid.ri
    t = data.t

    SigmaD = data.dust.Sigma
    SigmaG = data.gas.Sigma
    SigmaPlan = data.planetesimals.Sigma
    PlanMass = data.planetesimals.M
    eps = data.dust.eps

    cs = data.gas.cs
    delta = data.dust.delta.turb
    OmegaK = data.grid.OmegaK
    St = data.dust.St
    vK = OmegaK * r
    vFrag = data.dust.v.frag

    # Transformation of distribution
    a = np.mean(m[..., 1:] / m[..., :-1], axis=-1)
    dm = 2. * (a - 1.) / (a + 1.)
    sigmaD = SigmaD[..., :] / dm

    # Calculating limits

    # Fragmentation limit
    b = vFrag ** 2 / (delta * cs ** 2)
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings(
            'ignore',
            r'invalid value encountered in sqrt')
        St_fr = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b ** 2))

    # Drift limit
    p = SigmaG * OmegaK * cs

    _f = interp1d(np.log10(r), np.log10(p), fill_value='extrapolate')
    pi = 10. ** _f(np.log10(ri))
    gamma = np.abs(r / p * np.diff(pi) / np.diff(ri))
    St_dr = eps / gamma * (vK / cs) ** 2

    # Get limits
    sd_max = np.ceil(np.log10(sigmaD.max()))
    sg_max = np.ceil(np.log10(SigmaG.max()))

    # Text plotting
    center, width, frac, istart, iend = getRingStats(SigmaPlan, w)
    textstr = getText(PlanMassEarth, center, width, frac)
    titlestr = getTitle(z, w)

    fig, ax = plt.subplots()
    ax.set_title(titlestr, fontdict={'fontsize': fontsize})
    ax.text(0.04, 0.85, textstr, color="white", transform=ax.transAxes)


    pltcmap = ax.contourf(r / c.au, m, np.log10(sigmaD.T), levels=np.linspace(sd_max - 6, sd_max, 7), cmap="magma", extend="both")
    ax.contour(r / c.au, m, St.T, levels=[1.], colors="white", linewidths=2)
    ax.contour(r / c.au, m, (St - St_dr[..., None]).T, levels=[0.], colors="C2", linewidths=1.5)
    ax.contour(r / c.au, m, (St - St_fr[..., None]).T, levels=[0.], colors="C0", linewidths=1.5)
    ax.axhline(m[-1], color="#AAAAAA", lw=1, ls="--")
    ax.axvline(r[-1] / c.au, color="#AAAAAA", lw=1, ls="--")
    cbarcmap = plt.colorbar(pltcmap, ax=ax)
    cbarcmap.ax.set_ylabel("$\log\ \sigma$ [g/cm²]")
    ax.set_xlim(r[0] / c.au, 300)
    ax.set_ylim(m[0], m[-1])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Distance from star [au]")
    ax.set_ylabel("Particle mass [g]")
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(("%.0f"))
    ax.set_xticks([30, 60, 90, 150])

    fig.tight_layout()
    filename = outputDir + 'dist/d' + str(args.z)
    p = getPNG(filename)
    e = getPNG(filename)
    plt.savefig(p["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    plt.savefig(e["filename"], dpi=p["dpi"], bbox_inches=p["bbox"], pad_inches=p["pad"])
    # plt.savefig(filename + '.eps', format='eps')
    if args.show:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-s', action="store", dest="show", type=int, default=1, help="Show plot")
    arguments = parser.parse_args()
    main(arguments)
