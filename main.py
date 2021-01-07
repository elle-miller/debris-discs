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

try:
    from dustpy.simulation import Simulation
except:
    from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from dustpy.std.dust import MRN_distribution
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
    # Create instance
    s = Simulation()

    # Set initial conditions and initialise
    setInitConds(s, args, verbose=True)
    s.initialize()

    # Bind systoles and diastoles for planetesimal formation
    if args.planForm:
        setPlanForm(s)

    # Bind alpha bump function to create gas gap
    s.gas.alpha.updater.updater = alphaBumps
    s.update()

    # Bind initial gas profile and reinitialize
    s.gas.Sigma = initialGas(s, args.iniBumpPeakPos * c.au, args.amplitude, args.width, args.invertBump)
    s.dust.Sigma = MRN_distribution(s)
    s.update()

    # Specify where to put the data and simulation related info
    setSimulationParams(s, args)

    # Exlcude attributes to save space in simulation
    s.dust.backreaction.save = False
    s.dust.S.coag.save = False
    s.dust.p.frag.save = False
    s.dust.p.stick.save = False
    s.dust.v.rel.save = False
    if not args.dustEvolution:
        s.dust.H.save = False
        s.dust.SigmaFloor.save = False
        s.dust.v.rad.save = False
        s.dust.D.save = False

    # Run the simulation
    print("Evolving...")
    s.update()
    s.run()

    # Plot planel of results
    if args.panel:
        plot.ipanel(s.writer.datadir)


##################### OTHER FUNCTIONS #########################################

def setInitConds(s, args, verbose):
    """
    Set initial conditions of simulation model

    Input:
    ---------
    s: Instance of the Simulation class
    args: command line arguments
    verbose: bool
    """
    # Radial grid
    s.ini.grid.rmin = c.au * args.rmin
    s.ini.grid.rmax = c.au * args.rmax
    s.ini.grid.Nr = args.Nr
    # ri = np.logspace(np.log10(args.rmin), np.log10(args.rmax), args.Nr) * c.au
    # s.grid.ri = refinegrid(ri, (args.iniBumpPeakPos + 7 * args.width) * c.au)
    s.makegrids()

    # optional plotting of effect of refinement
    # fig, ax = plt.subplots()
    # ax.plot(s.grid.r / c.au, np.ones_like(s.grid.r), 'b*', label='r new')
    # ax.plot(s.grid.ri / c.au, np.ones_like(s.grid.ri), 'ro', label='ri refined')
    # ax.legend()
    # plt.show()

    # Mass grid
    # s.ini.grid.Nm = args.Nm
    # s.ini.grid.mmax = args.massMax  # default is 1e5

    # Gas
    alpha0 = args.alpha
    s.ini.gas.Mdisk = args.MdiskInMstar * c.M_sun  # args.MdiskInMstar set in solar masses so convert to grams
    s.ini.gas.alpha = alpha0 * np.ones_like(s.grid.r)

    # Dust (d2g ratio is dust.eps)
    # s.ini.dust.aIniMax = args.aIniMax
    s.ini.dust.vfrag = 1000.
    s.ini.dust.allowDriftingParticles=True

    # Star
    s.ini.star.M = c.M_sun * args.starmass

    # Bump
    bumpParams.init(A=args.amplitude, w=args.width, p=args.iniBumpPeakPos, v=args.bumpVelFactor, i=args.invertBump)

    if verbose:
        print("minyear = %d" % args.minyear)
        print("maxyear = %d" % args.maxyear)
        print("nsnap = %d" % args.nsnap)
        print("planForm = %s" % args.planForm)
        print("alpha = %f" % args.alpha)
        print("amplitude = %f" % args.amplitude)
        print("bumpVelFactor = %d" % args.bumpVelFactor)
        print("timeBumpForm = %f" % args.timeBumpForm)
        print("iniBumpPeakPos = %d" % args.iniBumpPeakPos)
        print("MdiskInMstar = %f" % args.MdiskInMstar)
        print("rmin = %d" % args.rmin)
        print("rmax = %d" % args.rmax)
        print("Nr = %d" % args.Nr)
        print("Nm = %d" % args.Nm)
        print("mmax = %d" % args.massMax)
        print("starmass = %f" % args.starmass)
        print("aIniMax = %f" % args.aIniMax)
        print("dustEvolution = %d" % args.dustEvolution)
        print("gasEvolution = %d" % args.gasEvolution)


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
    s.updater = ["star", "grid", "gas", "dust", "planetesimals"]
    s.planetesimals.updater = ["M"]

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

    # Simulation settings
    s.t.snapshots = np.logspace(args.minyear, args.maxyear, num=args.nsnap) * c.year
    print('snapshots=', len(s.t.snapshots))

    # Update the mixing params to match alpha (these used to be in "ini")
    s.dust.deltaRad = args.alpha  # radial particle diffusion
    s.dust.deltaTurb = args.alpha  # relative velocitiy turbulence
    s.dust.deltaVert = args.alpha  # vertical diffusion

    # Make the simulation shorter and smaller if not evolving gas or dust
    if not args.gasEvolution:
        del (s.integrator.instructions[1])

    if not args.dustEvolution:
        # Turn off all dust evolution
        s.dust.S.tot = 0.
        s.dust.S.updater = None
        s.dust.S.coag = 0
        s.dust.S.coag.updater = None
        s.dust.p.updater = None
        s.dust.v.rel.updater = None
        s.dust.v.frag.updater = None
        s.dust.kernel.updater = None
        s.dust.p.stick = 1.
        s.dust.p.frag = 0.
        s.dust.p.updater = None
        s.dust.Fi.adv = 0
        s.dust.Fi.adv.updater = None
        s.dust.v.rad = 0
        s.dust.v.rad.updater = None
        s.dust.Fi.diff = 0.
        s.dust.Fi.diff.updater = None
        s.dust.D.updater = None
        s.dust.S.hyd = 0.
        s.dust.S.hyd.updater = None
        s.dust.Fi.updater = None


class Bump:
    def __init__(self, args):
        self.alpha = args.alpha
        self.amplitude = args.amplitude
        self.position = args.iniBumpPeakPos
        self.width = args.width
        self.invert = args.invertBump


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
    parser.add_argument('-s', action="store", dest="minyear", type=float, default=3, help="Beginning year 10^x")
    parser.add_argument('-e', action="store", dest="maxyear", type=float, default=7, help="Ending year 10^x")
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
