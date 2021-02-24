from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from dustpy.std.dust import MRN_distribution
import numpy as np
from dustpy.std.gas import Hp
from matplotlib import pyplot as plt
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

# Using seaborn's style
fontsize = 14
labelsize = 14
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize, family='serif')
plt.rcParams['legend.fontsize'] = fontsize
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=fontsize)
fig, ax = plt.subplots()

# constants
au = 1.495978707e11
k = 1.380649e-23
mu = 2.3
mp = 1.67262192369e-27
G = 6.6743e-11
M = 1.988409870698051e30  # kg
R = 1391400000  # stellar radius (m)
T = 5772
phi = 0.05
B = (1.5 * k * T)/(mu * mp) * (R * np.sqrt(phi)/(G * M))**0.5
print("B=", B, " m/s")
print("or ", (B / au) * (1e6 * c.year), "au/Myr")

# Bump params
r_i = au * 90
alpha0 = 1e-3

time = 10e6 * c.year
f = np.linspace(0, 3, 1000)
alphas = [1e-3, 3e-4, 1e-4]
linetype = ['--', '-.', '-']
i = 0
a = 1e-3
ax.plot(f,  (r_i - a * B * time * f)/au, '--', label=r"$\alpha = 10^{-3}$")
a = 3e-4
ax.plot(f,  (r_i - a * B * time * f)/au, '-.', label=r"$\alpha = 3 \times 10^{-4}$")
a = 1e-4
ax.plot(f,  (r_i - a * B * time * f)/au, '-', label=r"$\alpha = 10^{-4}$")

# ax.set_title("Final bump distance after 10Myr, " + r"$\alpha$ = " + str(alpha0))
ax.set_xlabel("Velocity factor $f$")
ax.set_ylabel("Final distance [au]")
ax.grid(b=True)
# ax.plot(np.ones_like(f), np.linspace(0, r_i, 1000), 'r')
ax.set_ylim(0, r_i/au)
ax.set_xlim(0, 1)
filename = localDir + '/figplots/final_pos2'
ax.legend()
plt.savefig(filename+'.png', format='png', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.savefig(filename+'.eps', format='eps', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.show()
