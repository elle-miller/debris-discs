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

from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from functionsMovingBump import alphaBumps, initialGas
from functionsPlanFormation import setPlanetesimalFormation, dustSources, addPlanetesimals
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

    # Bind systoles and diastoles for planetesimal formation
    if args.planForm:
        setPlanForm(s)

    # Bind alpha bump function to create gas gap
    s.gas.alpha.updater.updater = alphaBumps
    s.update()

    # Bind initial gas profile and reinitialize
    #s.ini.gas.Sigma = initialGas(s, args.iniBumpPeakPos * c.au, args.amplitude, args.width)
    s.gas.Sigma = initialGas(s, args.iniBumpPeakPos * c.au, args.amplitude, args.width)
    s.ini.dust.allowDriftLimitedParticles = True
    s.initialize()

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
    s.run()

    # Plot planel of results
    if args.panel:
        plot.panel(s.writer.datadir)


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
    s.makegrids()

    # Mass grid
    s.ini.grid.Nm = args.Nm
    s.ini.grid.mmax = args.massMax  # default is 1e5
    if args.dustEvolution:
        # Emailed Stammler to see if necessary
        # Only enters for alpha0=1e-4
        if args.alpha < 3e-4:
            s.ini.grid.mmax = 1e14
        else:
            # This is where alpha0=1e-3 goes
            s.ini.grid.mmax = 1e8
    else:
        s.ini.grid.mmax = 1e1

    # Gas
    s.ini.gas.Mdisk = args.MdiskInMstar * c.M_sun  # args.MdiskInMstar set in solar masses so convert to grams
    s.ini.gas.alpha = args.alpha * np.ones_like(s.grid.r)

    # Dust (d2g ratio is dust.eps)
    s.ini.dust.aIniMax = args.aIniMax
    s.ini.dust.deltaRad = s.ini.gas.alpha  # radial particle diffusion
    s.ini.dust.deltaTurb = s.ini.gas.alpha  # relative velocitiy turbulence
    s.ini.dust.deltaVert = s.ini.gas.alpha  # vertical diffusion

    # Star
    s.ini.star.M = c.M_sun * args.starmass

    # Bump
    bumpParams.init(A=args.amplitude, w=args.width, p=args.iniBumpPeakPos, v=args.bumpVelFactor)

    s.initialize()

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

    # The functions we created above use fields that do not exist in a DustPy a-priori
    # We need to create the grids before the initialization
    s.dust.dSigmaDust = np.zeros((s.ini.grid.Nr, s.ini.grid.Nm))
    s.dust.SigmaPlan = np.ones(s.ini.grid.Nr) * 1.e-100

    # Calculate dust sources
    s.dust.updater.systole = setPlanetesimalFormation

    # Apply them
    s.dust.S.ext = dustSources(s)

    # Update after timestep
    s.dust.updater.diastole = addPlanetesimals



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
    s.t.snapshots = np.linspace(args.minyear, args.maxyear, num=args.nsnap, endpoint=True) * c.year
    print('snapshots=', len(s.t.snapshots))
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="outputDirNo", type=int, default=1, help="Simulation number")
    parser.add_argument('-l', action="store", dest="rmin", type=int, default=10, help="Inner radius limit (AU)")
    parser.add_argument('-c', action="store", dest="rmax", type=int, default=400, help="Outer radius limit (AU)")
    parser.add_argument('-s', action="store", dest="minyear", type=float, default=0, help="Beginning year")
    parser.add_argument('-e', action="store", dest="maxyear", type=float, default=10e6, help="Ending year")
    parser.add_argument('-n', action="store", dest="nsnap", type=int, default=100, help="Number of snapshots")
    parser.add_argument('-r', action="store", dest="Nr", type=int, default=100, help="Number of radial bins")
    parser.add_argument('-m', action="store", dest="Nm", type=int, default=120, help="Number of mass bins")
    parser.add_argument('-o', action="store", dest="starmass", type=float, default=1, help="Star Mass in solar units")
    parser.add_argument('-a', action="store", dest="alpha", type=float, default=0.001, help="Viscosity parameter")
    parser.add_argument('-b', action="store", dest="amplitude", type=float, default=10, help="log(Amplitude)")
    parser.add_argument('-w', action="store", dest="width", type=float, default=1., help="width factor")
    parser.add_argument('-t', action="store", dest="timeBumpForm", type=float, default=0, help="Time bump appears")
    parser.add_argument('-v', action="store", dest="bumpVelFactor", type=float, default=0, help="% of nominal")
    parser.add_argument('-p', action="store", dest="iniBumpPeakPos", type=int, default=90, help="Starting center (AU)")
    parser.add_argument('-i', action="store", dest="MdiskInMstar", type=float, default=0.1, help="Init disk mass (SM)")
    parser.add_argument('-g', action="store", dest="aIniMax", type=float, default=1e-4, help="Max initial dust size")
    parser.add_argument('-x', action="store", dest="massMax", type=float, default=1e5, help="Max mass")
    parser.add_argument('-1', action="store", dest="gasEvolution", type=int, default=1, help="Create bump via alpha")
    parser.add_argument('-2', action="store", dest="dustEvolution", type=int, default=1, help="Advect/diffus transport")
    parser.add_argument('-3', action="store", dest="panel", type=int, default=0, help="Plot panel output")
    parser.add_argument('-4', action="store", dest="planForm", type=int, default="1", help="Planetesimal formation")
    arguments = parser.parse_args()
    main(arguments)
