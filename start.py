from dustpy import Simulation
import dustpy.constants as c

import numpy as np

from simframe import Instruction
from simframe import schemes

import parameters as pars
import functions as f

import shutil
from os import path

"""
I think this is the main.py equivalent but if you want to do the refined radial grid stuff
"""


s = Simulation()
z = pars.z

# PARAMETERS

s.ini.dust.vfrag = pars.vfrag
s.ini.gas.alpha = pars.alpha
s.ini.gas.Mdisk = pars.Mdisk

# PLANET

# Group for planet with fields for location, velocity, width
s.addgroup("planet", description="Planetes quantities")
s.planet.addfield("r", pars.r_ini, description="Location [cm]")
s.planet.addfield("v", 0., description="Migration velocity [cm/s]")
# Approximate width (H/r = 0.05)
w = pars.f_w * 0.05 * s.planet.r
s.planet.addfield("w", w, description="Gap width [cm]")

# Updaters
s.planet.v.updater.updater = f.v_planet
s.planet.w.updater.updater = f.gap_width
# Update order
s.planet.updater = ["v", "w"]
# Differentiator
s.planet.r.differentiator = f.dr_planet
# Instruction
planet_instruction = Instruction(
    schemes.expl_1_euler, s.planet.r, description="Planet: explicit 1st-order Euler")

# RADIAL GRID

lmin = np.log10(pars.rmin)
lmax = np.log10(pars.rmax)
# We need a coarse grid and a fine grid.
# The coarse grid interfaces have to lie on the fine grid.
# Both are stored in fields.
ri_coarse = np.logspace(lmin, lmax, pars.Nr+1)
ri_fine = np.logspace(lmin, lmax, pars.f_ref*pars.Nr+1)
s.grid.addfield("ri_coarse", ri_coarse,
                description="Coarse grid interfaces [cm]")
s.grid.addfield("ri_fine", ri_fine, description="Fine grid interfaces [cm]")
# Approximate dust ring location for refinement
loc = s.planet.r + 2.5*s.planet.w
i = np.argmax(loc <= s.grid.ri_coarse) - 1
# This is the inde with the center of refinement in the coarse grid.
s.grid.addfield("i_ring", i, description="Index of dust ring in coarse grid")
s.grid.ri = f.get_grid(s, i)
# Initialize grids
s.makegrids()

# PLANETESIMALS

# Group for planetesimals with fields for mass and surface density
s.addgroup("planetesimals", description="Planetesimal quantities")
s.planetesimals.addfield("M", 0., description="Total planetesimal mass [g]")
s.planetesimals.addfield("Sigma", 1.e-100*np.ones(
    np.int(s.grid.Nr)), description="Planetesimal surface density [g/cm²]")
# Updater
s.planetesimals.M.updater.updater = f.M_plan
# Update order
s.planetesimals.updater = ["M"]
# Differentiator
s.planetesimals.Sigma.differentiator = f.dSigmaPlan
# Instruction
planetesimals_instruction = Instruction(schemes.expl_1_euler, s.planetesimals.Sigma,
                                        description="Planetesimals: explicit 1st-order Euler")

# UUpdate order of simulation frame
s.updater = ["star", "grid", "planet", "gas", "dust", "planetesimals"]


s.initialize()

# Setting snapshots
# 10 snapshots per decade until 1 Million years.
# Then 30 snapshots per decade.
s1 = 3.
s2 = 6.
N1 = np.int(10*(s2-s1))
s3 = 7.
N2 = np.int(30*(s3-s2))+1
snap1 = np.logspace(s1, s2, N1, endpoint=False)
snap2 = np.logspace(s2, s3, N2)
#s.t.snapshots = np.concatenate((snap1, snap2)) * c.year
s.t.snapshots = snap2 * c.year
print(s.t.snapshots / c.year )

# Output settings
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'
slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
s.writer.overwrite = True
if path.exists(localDir):
    outputDir = localDir + '/sims/' + str(z)
elif path.exists(slurmDir):
    outputDir = slurmDir + '/sims/' + str(z)
else:
    print("Output directory not found")
    exit(0)
s.writer.datadir = outputDir
if path.exists(outputDir):
    shutil.rmtree(outputDir)

# Adding integration instructions for planet location and planetesimal
# surface density to integrator
s.integrator.instructions.append(planet_instruction)
s.integrator.instructions.append(planetesimals_instruction)

# Loss term for dust
if pars.method == "tanh":
    s.dust.S.ext.updater.updater = f.S_ext_tanh
elif pars.method == "step":
    s.dust.S.ext.updater.updater = f.S_ext_step
else:
    raise(ValueError("Unknown method: {}".format(pars.method)))

# Set alpha updater
s.gas.alpha.updater.updater = f.alpha

# Update to get new alpha
s.update()

# Initializing deltas
s.dust.delta.rad = pars.delta_rad
s.dust.delta.turb = pars.delta_turb
s.dust.delta.vert = pars.delta_vert

# Set initial gas/dust densities with bump from new alpha
s.gas.Sigma *= pars.alpha / s.gas.alpha
s.dust.Sigma *= pars.alpha / s.gas.alpha[:, None]

# Final update before execution
s.update()

# Add systole to check if planet has left the simulation
# and for refining the radial grid.
s.updater.systole = f.systole

s.run()
