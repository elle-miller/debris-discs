import dustpy.constants as c

# file used for the alternate refined grid files: start.py, mainBasic.py

z = 318  # directory
Nr = 120

# PLANET
A_gap = 2.3  # Gap depth aplitude
f_v = 1.  # Migration velocity in units of viscous velocity
f_w = 1.  # Gap width in units of pressure scale height
r_ini = 90. * c.au  # Initial distance from star

# GRID

# The grid has a coarse grid and a fine grid. The coarse grid has Nr
# grid cells. The fine grid has f_ref*Nr grid cells.
# The actual grid on which the simulations runs is the coarse grid
# that is refined with 2N+1 fine grid cells centered around the dust ring.

# Number of grid cells on coarse grid

# Minimum and maximum grid extent
rmin = 10. * c.au
rmax = 400. * c.au
# Factor of refinement
f_ref = 3
# Range of refinement
N_ref = 7

# Disk parameters
Mdisk = 0.1 * c.M_sun

# Gas parameters
alpha = 1.e-3

# Dust parameters
vfrag = 1000.
delta_rad = 1.e-3
delta_turb = 1.e-3
delta_vert = 1.e-3

# Planetesimal formation
d2g_crit = 1.
zeta = 0.1
n = 0.03
method = "tanh"  # "tanh" or "step"
