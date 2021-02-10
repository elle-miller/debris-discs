from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from dustpy.std.dust import MRN_distribution
import numpy as np
from dustpy.std.gas import Hp
from matplotlib import pyplot as plt
au = 1.495978707e11
s = Simulation()

# Bump params
BumpPeakPos = au * 90
alpha0 = 1e-3
k = 1.38e-23
mu = 2.3
mp = 1.67e-27
G = 6.67e-11
M = 1.98e30
f = 0.1
time = 10e6
A = (mu * mp * np.sqrt(G * M)) / (160 * alpha0 * k * np.sqrt(au))
print("A = ", A)
f = np.linspace(0,3,1000)
x = 90 * np.exp(-time * f / A)
fig, ax = plt.subplots()
ax.plot(f, x)
plt.show()

# s.ini.dust.vfrag = 1000.
# s.ini.dust.allowDriftingParticles=True
#
# # Radial grid
# s.ini.grid.rmin = c.au * 10
# s.ini.grid.rmax = c.au * 400
# s.ini.grid.Nr = 200
# s.ini.grid.mmax = 1e8
# s.makegrids()
#
# s.initialize()
#
# iBumpPeakPos = (np.abs(s.grid.r - BumpPeakPos)).argmin()
#
# vr = -1.5 * BumpVelFactor * s.ini.gas.alpha * (s.gas.Hp[iBumpPeakPos] / s.grid.r[iBumpPeakPos]) ** 2 * \
#         s.grid.OmegaK[iBumpPeakPos] * s.grid.r[iBumpPeakPos]
# print(np.shape(vr))
# print(vr)
