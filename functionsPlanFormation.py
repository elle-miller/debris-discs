"""

Sets up rules for planetesimal formation. Mostly copied from DustPy except for the hyperbolic tangent code + boundaries.

Author: Elle Miller & Sebastian Marino (2021)
"""

import dustpy.std.sim as std_sim
import numpy as np
import time
import matplotlib.pyplot as plt


def M_plan(s):
    return (np.pi * (s.grid.ri[1:]**2 - s.grid.ri[:-1]**2) * s.planetesimals.Sigma[:]).sum()


def S_ext(s):

    # Planetesimal formation efficiency
    # zeta = 0.1

    # Midplane dust-to-gas ratio
    d2g_mid = s.dust.rho.sum(-1) / s.gas.rho

    # Old step function way
    # d2g_crit = 1.0
    # mask = np.where(d2g_mid >= d2g_crit, True, False)
    # ret = np.where(mask[:, None], -zeta * s.dust.Sigma * s.dust.St * s.grid.OmegaK[:, None], 0.)

    # New hyperbolic tangent way
    switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / s.bump.steep))
    ret = -s.bump.zeta * s.dust.Sigma * s.dust.St * s.grid.OmegaK[:, None] * switch[:, None]

    # Set to zero at boundaries - ret is Nr x Nm
    ret[0, :] = 0.
    ret[-10:, :] = 0.
    
    return ret


# The negative sum over all external source terms is the planetesimal formation rate
# We can write a source term function which we add as a derivative to the planetesiamal
# surface density for integration
def dSigmaPlan(s, x, Y):
    return -s.dust.S.ext.sum(-1)
