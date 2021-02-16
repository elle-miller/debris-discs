import numpy as np
from numpy import log, exp
from dustpy.std.gas import Hp


def initialGas(s):
    """
    Initial gas profile at t = 0. Must call after first initialize so all parameters are set
    """
    # Sigma without gap (tapered power law) and with gap
    r = s.grid.r
    SigmaGas_unperturbed = (r / s.ini.gas.SigmaRc) ** s.ini.gas.SigmaExp * np.exp(-(r / s.ini.gas.SigmaRc) ** (2+s.ini.gas.SigmaExp))
    SigmaGas_perturbed = SigmaGas_unperturbed * Gauss(s)

    # Normalize to get the right total disc mass
    M_gas = s.ini.gas.Mdisk / (s.ini.dust.d2gRatio + 1.)  # total gas mass in grams for given d2g
    normalization_factor = M_gas / np.trapz(2 * np.pi * r * SigmaGas_perturbed, x=r)
    iniGas = SigmaGas_perturbed*normalization_factor
    return iniGas


def alphaBumps(s):
    """
    Set turbulence values. This is the updater function for s.gas.alpha
    The current peak position is updated by calling BumpRadVel
    """
    # If bump pos < min exit the program, other return current peturbed alpha state
    s.bump.currentPeakPos = getPeakPosition(s)
    if s.bump.currentPeakPos < s.ini.grid.rmin:
        print("Exiting")
        exit()
    return s.ini.gas.alpha / Gauss(s)


def getPeakPosition(s):
    """
    Calculate the peak postion of the bump.
    """
    if s.t <= s.bump.timeStartMoving:
        peakPos = s.bump.iniPeakPos
    elif s.t > s.bump.timeStartMoving:
        peakPos = s.bump.currentPeakPos + BumpRadVel(s) * s.t.prevstepsize
    return peakPos


def BumpRadVel(s):
    """
    Radial velocity of a bump according to Armitage (2010, Eq. (7.32)).
    """
    r = s.grid.r
    iBumpPeakPos = (np.abs(r - s.bump.currentPeakPos)).argmin()
    return -1.5 * s.bump.f * s.ini.gas.alpha[0] * (s.gas.Hp[iBumpPeakPos] / s.grid.r[iBumpPeakPos]) ** 2 * \
        s.grid.OmegaK[iBumpPeakPos] * s.grid.r[iBumpPeakPos]


def Gauss(s):
    """
    Returns a Gaussian bump profile.
    """
    # Make sure bump at center of radial bin
    r = s.grid.r
    iBumpPeakPos = (np.abs(r - s.bump.currentPeakPos)).argmin()
    rpeak = r[iBumpPeakPos]
    w_gap = s.bump.width * Hp(s)  # Gas pressure scale height
    i = s.bump.invert
    return exp(-1 * i * log(s.bump.A) * exp(-0.5 * ((r - rpeak) / w_gap) ** 2))
