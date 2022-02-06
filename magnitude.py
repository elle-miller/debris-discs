import numpy as np
from matplotlib import pyplot as plt
from dustpy import constants as c
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

# This file prints the planetesimal surface density using the formula by Shibaike & ALibert (i think)

M_earth = 5.9722e24 * 1e3  # [g]

M_peb = 2e-5 * M_earth / c.year  # g/s
r = 50 * c.au  # cm
v_bump = 10 * c.au / c.year * 1e-6  # cm/s

sd_plan = M_peb / (2 * c.pi * r * v_bump)  # g s-1 cm-1 cm-1 s = g / cm2
print(sd_plan)
print("{t:.2e}".format(t=sd_plan))




