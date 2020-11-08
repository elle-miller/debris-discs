import dustpy.std.sim as std_sim
import numpy as np


##################### OTHER FUNCTIONS #########################################


# We transform some of the dust into planetesimals as soon as the dust to gas
# ratio is above a certain threshold
# This caculates the change rates for all dust species and the planetismal
# formation rate at each time step.
def setPlanetesimalFormation(s):
    # Just some helper variables to make it more readable
    OmegaK = s.grid.OmegaK  # Keplerian frequency          (Nr,)
    h = s.dust.H  # Dust scale heights           (Nr, Nm)
    SigmaDust = s.dust.Sigma  # Dust surface densities       (Nr, Nm)
    St = s.dust.St  # Dust Stokes numbers          (Nr, Nm)
    rhoGas = s.gas.rho  # Dust midplane volume density (Nr,)

    # First we calculate the total midplane volume density of the dust (Nr,)
    rhoDust = np.sum(SigmaDust / (np.sqrt(2. * np.pi) * h), axis=-1)  # <-- sum over mass axis
    # Then we calculate the midplane dust-to-gas ratio (Nr,)
    d2gMid = rhoDust / rhoGas

    # Now we need to check, where in the disk the critical midplane dust-to-gas is reached.
    d2gCrit = 1.  # This is the threshold we've chosen. Usually you would put this into a parameter file.
    trigger = np.where(d2gMid >= d2gCrit, True, False)  # (Nr,)

    # Here we calculate the change rates of the dust species. Everywhere where we trigger
    # streaming instability, we remove a fraction zeta of the dust per settling time scale.
    # The size-dependent settling time scale is given by 1/(St*OmegaK).
    zeta = 0.1  # Again hard-coded... Don't do this at home!
    dSigmaDust = np.where(trigger[:, None],  # Trigger expanded to (Nr, Nm)
                          -zeta * SigmaDust * St * OmegaK[:, None],  # The dust loss rate if triggered
                          0.)  # Else do nothing
    # For stability reasons at the boundaries we won't have planetesimal formation at the
    # first and last grid cell.
    dSigmaDust[0, :] = 0.
    dSigmaDust[-1, :] = 0.

    # The formation rate of planetesimals is then the negative sum over the mass dimension
    # of the dust loss rate
    dSigmaPlan = np.sum(-dSigmaDust, axis=-1)

    # Set the plan mass in the last 10 bins to zero
    dSigmaPlan[-10:] = 0.

    # We now store the change rates in our dust object to use them later.
    s.dust.dSigmaPlan = dSigmaPlan
    s.dust.dSigmaDust = dSigmaDust

    # For data analysis we also store the trigger, the midplane dust volume density, and
    # the midplane dust-to-gas ratio.
    s.dust.StrTrigger = trigger
    s.dust.rhoDust = rhoDust
    s.dust.d2gMid = d2gMid


# This function returns the source/loss terms we calculated before
def dustSources(s):
    return s.dust.dSigmaDust


def addPlanetesimals(s):
    #TODO: stepsize
    s.dust.SigmaPlan += s.t.stepsize * s.dust.dSigmaPlan