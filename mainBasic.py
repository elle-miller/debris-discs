#!/usr/bin/python3
# -*- coding: utf-8 -*-

###############################################################################
# main.py
# Author: Elle Miller
# Main file for program to evolve a protoplanetary disc under initial conditions
# To run:
# python main.py -flag1 val1 -flag2 val2 etc.
# Example script with dust evolution turned off
# python main.py -e 1e4 -n 5 -3 1 -2 0

# 82 just run plainly

try:
    from dustpy.simulation import Simulation
except:
    from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from functionsMovingBump import alphaBumps, initialGas
from functionsPlanFormation import M_plan, S_ext, dSigmaPlan
from simframe import Instruction, schemes
import matplotlib.pyplot as plt
import numpy as np
import argparse
import shutil
from os import path
import bumpParams


##################### MAIN FUNCTION ###########################################


def main(args):

    # Initial conditions
    alpha0 = 1.e-3
    vfrag = 1000.
    Rc = 60. * c.au
    Mdisk = 0.1 * c.M_sun
    r0 = 10. * c.au
    w = 2. * c.au
    A = 2.

    # Create instance and set necessary stuff
    s = Simulation()
    s.ini.gas.alpha = alpha0
    s.ini.gas.SigmaRc = Rc
    s.ini.gas.Mdisk = Mdisk
    s.ini.dust.vfrag = vfrag

    # Refine radial grid
    ri = np.logspace(0., 3., 100) * c.au
    s.grid.ri = refinegrid(ri, r0 + 3. * w)
    s.initialize()

    # Let simulation run until 1Myr to see PF
    s.t.snapshots = np.logspace(3., 6., 31) * c.year

    # alpha bump
    F = np.exp(-A * np.exp(-0.5 * (s.grid.r - r0) ** 2 / w ** 2))
    s.gas.alpha = alpha0 / F

    # Bind systoles and diastoles for planetesimal formation
    if args.planForm:
        setPlanForm(s)

    # Specify where to put the data and simulation related info
    setSimulationParams(s, args)

    # Run the simulation
    print("Evolving...")
    s.run()

    # Plot planel of results
    if args.panel:
        plot.ipanel(s.writer.datadir)


##################### OTHER FUNCTIONS #########################################


def setPlanForm(s):
    """
    Sets systole and diastole for planetesimal formation

    Input:
    ---------
    s: Instance of the Simulation class
    """

    # Adding planetesimals
    s.addgroup("planetesimals", description="Planetesimal quantities")
    s.planetesimals.addfield("M", 0., description="Total planetesimal mass [g]")
    s.planetesimals.addfield("Sigma", np.zeros_like(s.gas.Sigma),
                             description="Planetesimal surface density [g/cm²]")
    s.planetesimals.M.updater.updater = M_plan
    s.planetesimals.updater = ["M"]
    s.updater = ["star", "grid", "gas", "dust", "planetesimals"]

    # Adding planetesimal formation
    s.dust.S.ext.updater.updater = S_ext
    s.planetesimals.Sigma.differentiator = dSigmaPlan

    # Create an integration instruction and add it to the integrator
    inst_planetesimals = Instruction(
        schemes.expl_1_euler,
        s.planetesimals.Sigma,
        description="Planetesimals: explicit 1st-order Euler method")
    s.integrator.instructions.append(inst_planetesimals)
    s.update()


def setSimulationParams(s, args):
    """
    Sets output of simulation data and other simulation settings

    Input:
    ---------
    s: Instance of the Simulation class
    args: command line arguments
    """

    # Output settings
    localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'
    slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
    s.writer.overwrite = True
    if path.exists(localDir):
        outputDir = localDir + '/sims/' + str(args.outputDirNo)
    elif path.exists(slurmDir):
        outputDir = slurmDir + '/sims/' + str(args.outputDirNo)
    else:
        print("Output directory not found")
        return
    s.writer.datadir = outputDir
    if path.exists(outputDir):
        shutil.rmtree(outputDir)


def refinegrid(ri, r0, num=3):
    """Function to refine the radial grid

    Parameters
    ----------
    ri : array
        Radial grid
    r0 : float
        Radial location around which grid should be refined
    num : int, option, default : 3
        Number of refinement iterations

    Returns
    -------
    ri : array
        New refined radial grid"""
    if num == 0:
        return ri
    ind = np.argmin(r0 > ri) - 1
    indl = ind - num
    indr = ind + num + 1
    ril = ri[:indl]
    rir = ri[indr:]
    N = (2 * num + 1) * 2
    rim = np.empty(N)
    for i in range(0, N, 2):
        j = ind - num + np.int(i / 2)
        rim[i] = ri[j]
        rim[i + 1] = 0.5 * (ri[j] + ri[j + 1])
    ri = np.concatenate((ril, rim, rir))
    return refinegrid(ri, r0, num=num - 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="outputDirNo", type=int, default=1, help="Simulation number")
    parser.add_argument('-l', action="store", dest="rmin", type=int, default=10, help="Inner radius limit (AU)")
    parser.add_argument('-c', action="store", dest="rmax", type=int, default=400, help="Outer radius limit (AU)")
    parser.add_argument('-s', action="store", dest="minyear", type=float, default=0, help="Beginning year")
    parser.add_argument('-e', action="store", dest="maxyear", type=float, default=7, help="Ending year")
    parser.add_argument('-n', action="store", dest="nsnap", type=int, default=31, help="Number of snapshots")
    parser.add_argument('-r', action="store", dest="Nr", type=int, default=100, help="Number of radial bins")
    parser.add_argument('-m', action="store", dest="Nm", type=int, default=120, help="Number of mass bins")
    parser.add_argument('-o', action="store", dest="starmass", type=float, default=1, help="Star Mass in solar units")
    parser.add_argument('-a', action="store", dest="alpha", type=float, default=0.001, help="Viscosity parameter")
    parser.add_argument('-b', action="store", dest="amplitude", type=float, default=10, help="log(Amplitude)")
    parser.add_argument('-w', action="store", dest="width", type=float, default=1., help="width factor")
    parser.add_argument('-t', action="store", dest="timeBumpForm", type=float, default=0, help="Time bump appears")
    parser.add_argument('-v', action="store", dest="bumpVelFactor", type=float, default=0, help="% of nominal")
    parser.add_argument('-p', action="store", dest="iniBumpPeakPos", type=int, default=90, help="Starting center (AU)")
    parser.add_argument('-5', action="store", dest="MdiskInMstar", type=float, default=0.1, help="Init disk mass (SM)")
    parser.add_argument('-g', action="store", dest="aIniMax", type=float, default=1e-4, help="Max initial dust size")
    parser.add_argument('-x', action="store", dest="massMax", type=float, default=1e5, help="Max mass")
    parser.add_argument('-1', action="store", dest="gasEvolution", type=int, default=1, help="Create bump via alpha")
    parser.add_argument('-2', action="store", dest="dustEvolution", type=int, default=1, help="Advect/diffus transport")
    parser.add_argument('-3', action="store", dest="panel", type=int, default=0, help="Plot panel output")
    parser.add_argument('-4', action="store", dest="planForm", type=int, default="1", help="Planetesimal formation")
    parser.add_argument('-i', action="store", dest="invertBump", type=int, default="0", help="Invert the gauss bump")
    arguments = parser.parse_args()
    main(arguments)
