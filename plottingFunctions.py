from os import getcwd, path
from jobInfo import getJobParams
from dustpy import constants as c
import numpy as np

# Results are spead over two hard drives due to space
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]


def getDataDir(z):
    if path.exists(localDir + '/sims/' + str(z)):
        dataDir = localDir + '/sims/' + str(z)
    elif path.exists(localDirNew + '/' + str(z)):
        dataDir = localDirNew + '/' + str(z)
    else:
        print("Output directory not found")
        exit()
    return dataDir


def getTitle(z, writer):
    [alpha0, amplitude, position] = getJobParams(z)
    tMyrEnd = writer.read.sequence('t')[-1] / c.year * 1e-6

    # Wide scripts
    if z == 233:
        return r"$\alpha$: {a}, A: {A}, r$_p$: {p} au, f: 100\% at {t:.1f} Myr [{z}]".format(a=alpha0, A=amplitude, p=position,
                                                                                  t=tMyrEnd, z=str(z))
    elif z == 232 or z == 234:
        return r"$\alpha$: {a}, A: {A}, r$_p$: {p} au, f: 100\% after 1Myr, at {t:.1f} Myr [{z}]".format(a=alpha0, A=amplitude, p=position,
                                                                                  t=tMyrEnd, z=str(z))

    # Fragmentation scripts
    if 220 <= z <= 231:
        if np.mod(z, 2) == 0:
            titlestr = r"v$_f$: 1 m/s, $\alpha$: {a}, A: {A}, r$_p$: {p} au at {t:.1f} Myr [{z}]".format(a=alpha0, A=amplitude,
                                                                                       p=position, t=tMyrEnd, z=str(z))
        else:
            titlestr = r"v$_f$: 3 m/s, $\alpha$: {a}, A: {A}, r$_p$: {p} au at {t:.1f} Myr [{z}]".format(a=alpha0,
                                                                                                         A=amplitude,
                                                                                                         p=position,
                                                                                                         t=tMyrEnd,
                                                                                                         z=str(z))
        return titlestr
    stationary = False
    if z > 201:
        stationary = False
    elif z > 149 or z < 130:
        stationary = True
    if stationary:
        titlestr = r"$\alpha$: {a}, A: {A}, r$_p$: {p} au at {t:.1f} Myr [{z}]".format(a=alpha0, A=amplitude,
                                                                                       p=position, t=tMyrEnd, z=str(z))
    else:
        titlestr = r"$\alpha$: {a}, A: {A}, f: {p}\% at {t:.1f} Myr [{z}]".format(a=alpha0, A=amplitude, p=position,
                                                                                  t=tMyrEnd, z=str(z))
    return titlestr


def getText(p, c1, w1, f1):
    ptot = f"{p:.0f}"
    textstr = r"$\mathcal{M}$: " + str(ptot) + r" M$_{\oplus}$" + "\n"
    center = f"{c1:.0f}"
    width = f"{w1:.0f}"
    frac = f"{f1:.1f}"
    if w1 > 0.05:
        textstr += "c: " + center + " au, w: " + width + " au, f: " + frac
    return textstr


def getRingStats(SigmaPlan, writer, verbose=False):
    """
    :return: center, width & fractional width of a given planetesimal surface density profile
    """
    rInt = writer.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]

    # Reset index positions
    istartRingP = 0
    iendRingP = 0
    index = 0
    minVal = 1e-3

    # Loop through each radial bin, locating index positions of start and end ring
    for k in SigmaPlan:
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
        # Left interface
        startRing = rInt[-1, istartRingP] / c.au
        # Right hand side of bin, so +1
        endRing = rInt[-1, iendRingP + 1] / c.au

    c1 = (endRing + startRing) / 2
    w1 = endRing - startRing
    f1 = w1 / c1
    if verbose:
        print("Ring Start: %.1f" % startRing)
        print("Ring End: %.1f" % endRing)
        print("Ring Center: %.1f AU" % c1)
        print("Ring Width: %.1f AU" % w1)
        print("Ring Fractional Width: %.2f" % f1)
    return c1, w1, f1, istartRingP, iendRingP


def getRingMass(RingDustTot, istartRingP, iendRingP, writer):
    # If there are no rings, no ring dust
    if istartRingP == 0:
        RingDustTot[:] = 0

    # Else, mask dust not within the indices
    else:
        Nr = np.size(RingDustTot)
        for k in range(Nr):
            if k not in range(istartRingP, iendRingP + 1):
                RingDustTot[k] = 0

    rInt = writer.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]
    RingDiscMass = np.sum(np.pi * (rInt[-1, 1:] ** 2. - rInt[-1, :-1] ** 2.) * RingDustTot[:], axis=0) / M_earth
    return RingDiscMass

def getEPS(filename):
    dict = {
        "filename": filename + ".eps",
        "format": "eps",
        "dpi": 300,
        "bbox": "tight",
        "pad": 0.1
    }
    return dict

def getPNG(filename):
    dict = {
        "filename": filename + ".png",
        "format": "png",
        "dpi": 300,
        "bbox": "tight",
        "pad": 0.1
    }
    return dict
