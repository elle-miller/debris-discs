import numpy as np
from numpy import log, exp
from dustpy.std.gas import Hp
from dustpy import constants as c
from matplotlib import pyplot as plt
from dustpy.std.sim import dt
import bumpParams


def initialGas(s, IniBumpPeakPos, A, width):
    """
    Initial gas profile at t = 0. Must call after first initialize so all parameters are set

    Input:
    ---------
    IniBumpPeakPos: float
    position of bump peak
    A: float
    amplitude
    width: float
    bump width scale factor
    """
    r = s.grid.r
    SigmaGas = (r / s.ini.gas.SigmaRc) ** -s.ini.gas.SigmaExp * np.exp(-(r / s.ini.gas.SigmaRc) ** (2-s.ini.gas.SigmaExp))
    M_gas = s.ini.gas.Mdisk * s.ini.star.M / (s.ini.dust.d2gRatio + 1.)  # total gas mass for given d2g
    s.ini.gas.Sigma0 = M_gas / np.trapz(2 * np.pi * r * SigmaGas, x=r)
    BumpPeakPos = getPeakPosition(s, IniBumpPeakPos=IniBumpPeakPos, TimeBumpForm=0, BumpVelFactor=0)
    iniGas = Gauss(s, r, BumpPeakPos, A, width) * s.ini.gas.Sigma0 * SigmaGas
    return iniGas


def Gauss(s, r, BumpPeakPos, A, width):
    """
    Gaussian bump.

    Input:
    ---------
    A: float
    amplitude
    r: array
    radial grid
    rpeak: float
    peak position of bump
    width: float
    bump width scale factor
    """
    # Make sure bump at center of radial bin
    iBumpPeakPos = (np.abs(r - BumpPeakPos)).argmin()
    rpeak = r[iBumpPeakPos]
    w_gap = width * Hp(s)  # Gas pressure scale height
    return exp(-log(A) * exp(-0.5 * ((r - rpeak) / w_gap) ** 2))


def BumpRadVel(s, iBumpPeakPos, BumpVelFactor):
    """
    Radial velocity of a bump according to Armitage (2010, Eq. (7.32)).
    Input:
    --------
    iBumpPeakPos: integer
        bin closest to current gas bump peak position
    """
    # print(BumpVelFactor)
    # print(repr(s.ini.gas.alpha))
    # print(s.gas.Hp[iBumpPeakPos])
    # print(s.grid.r[iBumpPeakPos])
    # print(s.grid.OmegaK[iBumpPeakPos])
    return -1.5 * BumpVelFactor * s.ini.gas.alpha[0] * (s.gas.Hp[iBumpPeakPos] / s.grid.r[iBumpPeakPos]) ** 2 * \
        s.grid.OmegaK[iBumpPeakPos] * s.grid.r[iBumpPeakPos]


def getPeakPosition(s, IniBumpPeakPos, TimeBumpForm, BumpVelFactor):
    """
    Calculate the peak postion of the bump.

    Input:
    --------
    IniBumpPeakPos: float
        initial peak postion of bump
    vRadBump: float
        radial velocity of the moving gas bump at the peak position
    TimeBumpForm: float
        time at which a gas bump is forced
    """
    # Define and save the peak position
    if s.t <= TimeBumpForm:
        s.gas.BumpPeakPos = IniBumpPeakPos
    elif s.t > TimeBumpForm:
        iBumpPeakPos = (np.abs(s.grid.r - s.gas.BumpPeakPos)).argmin()
        s.gas.BumpPeakPos = IniBumpPeakPos
        s.gas.BumpPeakPos += BumpRadVel(s, iBumpPeakPos=iBumpPeakPos, BumpVelFactor=BumpVelFactor) * s.t.prevstepsize

        # print("iBumpPeakPos= ", iBumpPeakPos)
        # print("BumpPeakPos= ", IniBumpPeakPos/c.au)
        # print("Previous step size= ",repr(s.t.prevstepsize))
        # print("BumpRadVel= ", repr(BumpRadVel(s, iBumpPeakPos=iBumpPeakPos, BumpVelFactor=BumpVelFactor)))
        # print("New BumpPeakPos= ", repr(s.gas.BumpPeakPos/c.au))

    return s.gas.BumpPeakPos


def renormalizeGasProfile(s, M_gas, IniBumpPeakPos, BumpVelFactor, A, width, TimeBumpForm):
    """
    For a non-evolving disk, the gas density is normalized in every timestep to guarantee mass conservation.
    Input:
    --------
    M_gas: float
        total gas disk mass [g]
    width: float
    bump width scale factor
    """

    if s.t >= TimeBumpForm:
        r = s.grid.r
        s.gas.Sigma = (r / s.ini.gas.SigmaR0) ** -s.ini.gas.SigmaExp * np.exp(
            -(r / s.ini.gas.SigmaR0) ** (2 - s.ini.gas.SigmaExp))
        BumpPeakPos = getPeakPosition(s, IniBumpPeakPos, TimeBumpForm, BumpVelFactor)
        s.gas.Sigma += Gauss(s, r, BumpPeakPos, A, width)
        SigmaNormConst = M_gas / np.trapz(2 * np.pi * s.grid.r * s.gas.Sigma, x=s.grid.r)
        s.gas.Sigma = SigmaNormConst * s.gas.Sigma
    return s.gas.Sigma


def alphaBumps(s, IniBumpPeakPos, A, width, TimeBumpForm, BumpCreatedViaAlpha, BumpVelFactor):
    """
    Set turbulence values. This function is called in the code and
    will be used to call the other funtions that change the gas density.

    Input:
    --------
    IniBumpPeakPos: float
    initial position of bump
    A: float
    amplitude of bump
    TimeBumpForm: float
    time after which a bump is allowed to form
    width: float
    bump width scale factor
    """
    r = s.grid.r
    bumpyAlpha = s.ini.gas.alpha * np.ones_like(r)
    M_gas = s.ini.gas.Mdisk * s.ini.star.M / (s.ini.dust.d2gRatio + 1.)  # total gas mass for given d2g

    if not BumpCreatedViaAlpha:
        s.gas.Sigma = renormalizeGasProfile(s, M_gas, IniBumpPeakPos, BumpVelFactor, A, width, TimeBumpForm)
    else:
        # If the peak position of the bump goes beyond the set minimum radius, exit the program
        BumpPeakPos = getPeakPosition(s, IniBumpPeakPos, TimeBumpForm, BumpVelFactor)
        if BumpPeakPos < s.ini.grid.rmin:
            print("Exiting")
            exit()
        else:
            bumpyAlpha = s.ini.gas.alpha / Gauss(s, r, BumpPeakPos, A, width)

    # Show a plot of alpha wrt radius
    fig, ax = plt.subplots()
    ax.loglog(r/c.au, bumpyAlpha, label="Alpha")
    plt.show()

    return bumpyAlpha


# Alpha bump with one input argumen
def alphaBumps2(s):
    """
    Set turbulence values. This function is called in the code and
    will be used to call the other funtions that change the gas density.

    Input:
    --------
    IniBumpPeakPos: float
    initial position of bump
    A: float
    amplitude of bump
    TimeBumpForm: float
    time after which a bump is allowed to form
    width: float
    bump width scale factor
    """
    IniBumpPeakPos = bumpParams.position*c.au
    A = bumpParams.amplitude
    width = bumpParams.width
    BumpVelFactor = bumpParams.velocity
    BumpCreatedViaAlpha = True
    TimeBumpForm = 0
    r = s.grid.r
    bumpyAlpha = s.ini.gas.alpha * np.ones_like(r)
    M_gas = s.ini.gas.Mdisk * s.ini.star.M / (s.ini.dust.d2gRatio + 1.)  # total gas mass for given d2g

    if not BumpCreatedViaAlpha:
        s.gas.Sigma = renormalizeGasProfile(s, M_gas, IniBumpPeakPos, BumpVelFactor, A, width, TimeBumpForm)
    else:
        # If the peak position of the bump goes beyond the set minimum radius, exit the program
        BumpPeakPos = getPeakPosition(s, IniBumpPeakPos, TimeBumpForm, BumpVelFactor)
        if BumpPeakPos < s.ini.grid.rmin:
            print("Exiting")
            exit()
        else:
            bumpyAlpha = s.ini.gas.alpha / Gauss(s, r, BumpPeakPos, A, width)
    return bumpyAlpha
