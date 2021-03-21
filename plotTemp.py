from dustpy import plot
from dustpy import Simulation
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import numpy as np
import argparse
import matplotlib.pyplot as plt
from movieBump import movieBump
from dustpy.plot import panel
from os import path, getcwd
from matplotlib.ticker import ScalarFormatter
from plottingFunctions import *
from dustpy import hdf5writer as w


s = Simulation()
s.ini.grid.rmin = c.au * 10
s.ini.grid.rmax = c.au * 400
s.ini.grid.Nr = 200
s.ini.grid.mmax = 1e8
s.makegrids()
s.ini.gas.Mdisk = 0.1 * c.M_sun
s.ini.gas.alpha = 0.001 * np.ones_like(s.grid.r)
s.initialize()

tstar = s.star.T
rstar = s.star.R

phi = 0.05
#print(tstar, rstar* 1e-5, phi)
temp = tstar * (rstar * phi ** 0.5) ** 0.5 / np.sqrt(c.au)
print("{t:.1f}".format(t=temp))
exit(0)

# Plot temperature profile
fig, ax = plt.subplots()

ax.plot(s.grid.r / c.au, s.gas.T, label="Real profile")
ax.plot(s.grid.r / c.au, 263 * (c.au/s.grid.r)**0.5, label="Approx ~263")
ax.plot(s.grid.r / c.au, 260 * (c.au/s.grid.r)**0.5, label="Approx ~260")
ax.legend()
ax.set_title("Temperature profile of disc")
ax.set_ylabel("Temperature [K]")
ax.set_xlabel("Radius [au]")
plt.show()