import numpy as np
from scipy.interpolate import interp1d
import sys

import parameters as pars
import dustpy.constants as c


def gap_width(s):
    """Function returns the gap width. Factor of pressure scale height
    at planet location."""
    i = (np.abs(s.grid.r - s.planet.r)).argmin()
    return pars.f_w * s.gas.Hp[i]


def v_planet(s):
    """Function returns migration velocity of planet."""
    r = s.grid.r
    i = (np.abs(r - s.planet.r)).argmin()
    return -1.5 * pars.f_v * pars.alpha * (s.gas.Hp[i] / r[i])**2 * s.grid.OmegaK[i] * r[i]


def dr_planet(s, x, Y):
    """Differential equation of planet location.
    Simply the migration speed."""
    return s.planet.v


def alpha(s):
    """Modified alpha profile to create gap."""
    r = s.grid.r
    i = (np.abs(r - s.planet.r)).argmin()
    return pars.alpha / np.exp(-pars.A_gap * np.exp(-0.5*(r-r[i])**2/s.planet.w**2))


def M_plan(s):
    """Function returns the total planetesimal mass."""
    return (s.grid.A*s.planetesimals.Sigma).sum(-1)


def S_ext_step(s):
    """Function returns the loss terms of the dust surface density.
    Step function that turns on when the dust-to-gas ratio is larger than one."""

    # Midplane dust-to-gas ratio
    d2g_mid = s.dust.rho.sum(-1) / s.gas.rho

    # Mask that defines if planetesimal formation is triggered
    mask = np.where(d2g_mid >= pars.d2g_crit, True, False)

    # Change in dust surface densities
    ret = np.where(mask[:, None], -pars.zeta*s.dust.Sigma *
                   s.dust.St * s.grid.OmegaK[:, None], 0.)
    # Set to zero at boundaries
    ret[0, :] = 0.
    ret[-1, :] = 0.

    return ret


def S_ext_tanh(s):
    """Function returns the loss terms of the dust surface density.
    Smooth transition in dust-to-gas ratio"""

    # Midplane dust-to-gas ratio
    d2g_mid = s.dust.rho.sum(-1) / s.gas.rho

    # New hyperbolic tangent way
    switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / pars.n))
    ret = -pars.zeta * s.dust.Sigma * s.dust.St * \
        s.grid.OmegaK[:, None] * switch[:, None]

    # Set to zero at boundaries
    ret[0, :] = 0.
    ret[-1, :] = 0.

    return ret


def dSigmaPlan(s, x, Y):
    """Differential equation of planetesimal surface density.
    Simply the negative sum of all dust loss terms."""
    return -s.dust.S.ext.sum(-1)


def systole(s):
    """Systole that is executed once per time step.
    The function checks if the planet is still within the grid and it
    will refine the grid if necessary."""

    # Check if planet is close to inner edge
    planetexit(s)

    # Adjust grid if planet moved
    adjustgrid(s)


def planetexit(s):
    """Function terminates simulation if planet is close to inner edge.
    The function writes an output with the last simulation state,
    even if no snapshot should be written."""
    if s.planet.r-2.5*s.planet.w <= s.grid.ri[0]:
        # Save snapshot
        i = np.argmax(s.t <= s.t.snapshots) + 1
        s.writeoutput(i)
        sys.exit("The planet has left the building.")


def adjustgrid(s):
    """This function adjusts the grid if the dust ring location
    moved changed the cell in the coarse grid."""

    # Get new dust ring location
    loc = s.planet.r + 2.5*s.planet.w
    # Index of ring location and new index of ring location
    i = s.grid.i_ring
    i_n = np.argmax(loc <= s.grid.ri_coarse) - 1

    # Do nothing if index has not changed (ie. the ring did not change grid cell)
    if i == i_n:
        return

    # If the index changed, set new grid
    set_new_grid(s, i_n)


def set_new_grid(s, i_n):
    """Function creates new grid and calculates the quantities on the grid
    such that mass is conserved."""

    # Storing required quantities in short variables.
    i = np.int(s.grid.i_ring)
    i_d = i_n - i  # Difference in ring index
    r = s.grid.r
    ri = s.grid.ri
    ric = s.grid.ri_coarse
    rif = s.grid.ri_fine
    ri_n = get_grid(s, i_n)  # New gridcell interfaces
    r_n = 0.5 * (ri_n[1:] + ri_n[:-1])
    A_n = np.pi*(ri_n[1:]**2 - ri_n[:-1]**2)
    m = pars.f_ref
    N = pars.N_ref

    # Calculate quantities on new grid.
    s.planetesimals.Sigma[...] = get_quantity(
        i, i_d, m, N, r, r_n, ri, ri_n, ric, rif, s.planetesimals.Sigma)
    s.gas.Sigma[...] = get_quantity(
        i, i_d, m, N, r, r_n, ri, ri_n, ric, rif, s.gas.Sigma)
    # That would be better without loop.
    for k in range(np.int(s.grid.Nm)):
        s.dust.Sigma[..., k] = get_quantity(
            i, i_d, m, N, r, r_n, ri, ri_n, ric, rif, s.dust.Sigma[..., k])

    # Store new grid quantities.
    s.grid.r[...] = r_n
    s.grid.ri[...] = ri_n
    s.grid.A[...] = A_n
    s.grid.i_ring = i_n


def get_index(ic, ria, ric):
    """Function returns index that an index on the coarse grid
    has on the actual adaptive grid.

    Parameters
    ----------
    ic : integer
        Index on coarse grid
    ria : array
        Adaptive grid
    ric : array
        Coarse grid

    Returns
    -------
    ia : integer
        Index on adaptive grid"""
    return np.isclose(ria, ric[ic]).argmax()


def to_coarse(rif, arr):
    """Function transforms a quantity on the fine grid onto the coarse grid,
    such that mass is conserved.

    Parameters
    ----------
    rif : array(N+1)
        Part of the fine grid cell interfaces
    arr : array(N)
        Quantity on the fine grid cell centers

    Returns
    -------
    arr_c : float
        Quantity on the coarse grid"""
    return ((rif[1:]**2-rif[:-1]**2)*arr[...]).sum(-1) / (rif[-1]**2-rif[0]**2)


def to_fine(r_before, sig_before, r_after, ri_after):
    """Function transforms quantity from coarse grid onto fine grid with
    interpolation and normalization, such that mass is conserved.

    Parameters
    ----------
    r_before : array(3)
        Values of coarse grid cell centers including left and right cell
    sig_before: array(3)
        Values of quantity at r_before
    r_after : array(N)
        Grid cell centers of new fine grid
    ri_after : array(N)
        Grid cell interfaces of new fine grid

    Returns
    -------
    sig_after : array(N)
        New values of quantity on fine grid"""

    # Total mass in grid cell
    M_before = (ri_after[-1]**2 - ri_after[0]**2)*sig_before[1]
    # Interpolation onto new grid
    f = interp1d(r_before, sig_before*r_before)
    sig_after = f(r_after)/r_after
    # Mass after the interpolation
    M_after = ((ri_after[1:]**2 - ri_after[:-1]**2) * sig_after).sum()
    # Normalization to conserve mass
    sig_after *= M_before/M_after
    return sig_after


def get_quantity(i, i_d, m, N, r, r_n, ri, ri_n, ric, rif, sig):
    """Function returns quantity on the new grid.

    Parameters
    ----------
    i : integer
        Old index of dust ring on coarse grid
    i_d : integer
        Change in index
    m : integer
        Refinement factor
    N : integer
        Refinement range
    r : array(Nr)
        Old radial grid cell centers
    r_n : array(Nr)
        New radial grid cell centers
    ri : array(Nr+1)
        Old radial grid cell interfaces
    ri_n : array(Nr+1)
        New radial grid cell interfaces
    ric : array
        Coarse grid cell interfaces
    rif : array
        Fine grid cell interfaces
    sig : array(Nr)
        Qantity on old grid cell centers

    Returns
    -------
    sig_n : array(Nr)
        Quantity on new grid"""
    # Initialize with zeros
    sig_n = np.zeros_like(sig)

    # Coarse index movement to left (zero or negative)
    i_neg = np.minimum(i_d, 0)
    # Coarse index movement to right (zero or positive)
    i_pos = np.maximum(0, i_d)

    # The inner and outer coarse region that remains unchanged

    # Last inner unchanged coarse index
    k_c_in = i - N + i_neg
    # First outer unchanged coarse index
    k_c_out = i + N + 1 + i_pos
    # Convert to adaptive grid indices
    k_in = k_c_in  # The inner indices are always indentical, because they are before the refinement
    k_in_n = k_c_in
    k_out = get_index(k_c_out, ri, ric)
    k_out_n = get_index(k_c_out, ri_n, ric)
    # Simply copy the unchanged coarse values
    sig_n[:k_in_n] = sig[:k_in]
    sig_n[k_out_n:] = sig[k_out:]

    # The central fine region that remains unchanged

    # First coarse index of unchanged fine grid
    k_f_in = i-N+i_pos
    # Last coarse index of unchanged fine grid
    k_f_out = i+N+i_neg
    # Convert to adaptive grid indices
    k_in = get_index(k_f_in, ri, ric)
    k_out = get_index(k_f_out+1, ri, ric)
    k_in_n = get_index(k_f_in, ri_n, ric)
    k_out_n = get_index(k_f_out+1, ri_n, ric)
    # Copy unchanged fine data
    sig_n[k_in_n:k_out_n] = sig[k_in:k_out]

    # The fine region that changed

    for i in range(np.abs(i_d)):
        k_in = get_index(k_c_in+i, ri, ric)
        k_in_n = get_index(k_c_in+i, ri_n, ric)
        k_out = get_index(k_c_out-1-i, ri, ric)
        k_out_n = get_index(k_c_out-1-i, ri_n, ric)
        if i_d > 0:  # Grid moved to the right
            # Inner region
            sig_n[k_in_n] = to_coarse(ri[k_in:k_in+m+1], sig[k_in:k_in+m])
            # Outer region
            r_b = r[k_out-1:k_out+2]
            sig_before = sig[k_out-1:k_out+2]
            r_a = r_n[k_out_n:k_out_n+m]
            ri_a = ri_n[k_out_n:k_out_n+m+1]
            sig_n[k_out_n:k_out_n+m] = to_fine(r_b, sig_before, r_a, ri_a)
        else:  # Grid moved to the left
            # Inner region
            r_b = r[k_in-1:k_in+2]
            sig_before = sig[k_in-1:k_in+2]
            r_a = r_n[k_in_n:k_in_n+m]
            ri_a = ri_n[k_in_n:k_in_n+m+1]
            sig_n[k_in_n:k_in_n+m] = to_fine(r_b, sig_before, r_a, ri_a)
            # Outer region
            sig_n[k_out_n] = to_coarse(
                ri[k_out:k_out+m+1], sig[k_out:k_out+m])

    return sig_n


def get_grid(s, i):
    """Function returns grid with refinement centered on i."""
    # Helpers
    N = pars.N_ref
    m = pars.f_ref
    ri_c = s.grid.ri_coarse
    ri_f = s.grid.ri_fine

    # Inner and outer coarse grid and middle fine grid
    ri_in = ri_c[:i-N]
    ri_out = ri_c[i+N+1:]
    ri_mid = ri_f[(i-N)*m:(i+N+1)*m]

    return np.concatenate((ri_in, ri_mid, ri_out))
