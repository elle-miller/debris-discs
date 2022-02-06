from dustpy import constants as c
import numpy as np

# This file calculates the value of B from the paper

# constants in SI
au = 1.495978707e11
k = 1.380649e-23
mu = 2.3  # kg
mp = 1.67262192369e-27
G = 6.6743e-11
M = 1.988409870698051e30  # kg
R = 1391400000  # stellar radius (m)
T = 5772
phi = 0.05
B = (1.5 * k * T)/(mu * mp) * (R * np.sqrt(phi)/(G * M))**0.5
print("B=", B, " m/s")
B = B / au * c.year
print("or {B:.2e} au/Myr".format(B=B * 1e6))
B = 10 ** 4 / 1e6
