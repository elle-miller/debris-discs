import dustpy.std.sim as std_sim
import numpy as np


##################### OTHER FUNCTIONS #########################################

def M_plan(s):
    return (np.pi * (s.grid.ri[1:]**2 - s.grid.ri[:-1]**2) * s.planetesimals.Sigma[:]).sum()


def S_ext(s):

    # Critical dust-to-gas ratio
    d2g_crit_max = 1.
    d2g_crit_cent = 0.75
    d2g_crit_min = 0.5

    # Planetesimal formation efficiency
    zeta = 0.1

    # Midplane dust-to-gas ratio
    d2g_mid = s.dust.rho.sum(-1) / s.gas.rho

    # Mask that defines if planetesimal formation is triggered
    # mask = np.where(d2g_mid >= d2g_crit, True, False)

    # Change in dust surface densities
    # ret = np.where(mask[:, None], -zeta*s.dust.Sigma * s.dust.St * s.grid.OmegaK[:, None], 0.)
    switch = 0.5 * (1 + np.tanh((d2g_mid - 0.75) / 0.25))
    # ret = -zeta*s.dust.Sigma * s.dust.St * s.grid.OmegaK[:, None] * 0.5 * (1 + np.tanh((d2g_mid-0.75)/0.25))
    ret = -zeta * s.dust.Sigma * s.dust.St * s.grid.OmegaK[:, None] * switch[:, None]

    # Set to zero at boundaries
    ret[0, :] = 0.
    ret[-5, :] = 0.

    return ret


# The negative sum over all external source terms is the planetesimal formation rate
# We can write a source term function which we add as a derivative to the planetesiamal
# surface density for integration
def dSigmaPlan(s, x, Y):
    return -s.dust.S.ext.sum(-1)
