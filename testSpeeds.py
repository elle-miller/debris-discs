from dustpy.simulation import Simulation
from dustpy import constants as c
from dustpy import plot
from dustpy.std.dust import MRN_distribution
import numpy as np
from dustpy.std.gas import Hp
from matplotlib import pyplot as plt
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

# Using seaborn's style
width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
labelsize = 14
plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches
fig, ax = plt.subplots()

colors = ["darkmagenta", "orchid"]

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
B = B / au * c.year

time = np.linspace(0, 10e6, 1000)
f = 1
a = 1e-3
ax.plot(time/1e6,  (a * B * time * f), '-', color=colors[0], linewidth=2, label=r"$\alpha = 10^{-3}$")
a = 1e-4
ax.plot(time/1e6,  (a * B * time * f), '-', color=colors[1], linewidth=2, label=r"$\alpha = 10^{-4}$")

# ax.set_title("Final bump distance after 10Myr, " + r"$\alpha$ = " + str(alpha0))
ax.set_xlabel("Time [Myr]")
ax.set_ylabel("Distance travelled [au]")
ax.grid(b=True)
ax.set_ylim(0, 105)
ax.set_xlim(0, 10)
filename = localDir + '/figplots/width_vs_time'
ax.legend(framealpha=1.0)
plt.savefig(filename+'.png', format='png', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.savefig(filename+'.eps', format='eps', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.show()

yeet = True
endtime = 7.3564225445964215
if yeet:
    fig, ax = plt.subplots()
    time = endtime * 1e6
    f = np.linspace(0, 100, 10000)
    a = 1e-3
    ax.plot(f, 90 - (a * B * time * f/100), '-', color=colors[0], linewidth=2, label=r"$\alpha = 10^{-3}$")
    a = 1e-4
    ax.plot(f, 90 - (a * B * time * f/100), '-', color=colors[1], linewidth=2, label=r"$\alpha = 10^{-4}$")
    ax.vlines(10, 10, 100, ls='--', color='pink')
    ax.vlines(30, 10, 100, ls='--', color='pink')
    ax.vlines(100, 10, 100, ls='--', color='pink')
    # ax.set_title("Final bump position after 7.4 Myr")
    ax.legend(framealpha=1.0)
    ax.set_xlabel(r"Bump velocity $f$ as $\%$ of nominal")
    ax.set_ylabel("Peak position [au]")
    ax.set_title("Theoretical bump position after 7.4 Myr", fontdict={'fontsize': fontsize})
    ax.grid(b=True)
    ax.set_ylim(10, 100)
    ax.set_xlim(0, 100)
    plt.savefig(filename + 'v2.png', bbox_inches='tight', pad_inches=0.05, dpi=300)
    plt.show()

