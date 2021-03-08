import numpy as np
from matplotlib import pyplot as plt
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

# Figure stuff
width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14
labelsize = 14

plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches
fig, ax = plt.subplots()

rgb = np.array([(158, 201, 226), (60, 147, 194), (13, 74, 112)]) / 255  # Mono-hue blue
colors = np.array([rgb[2], rgb[1], rgb[0]])
colors = np.array(["midnightblue", "royalblue", "lightsteelblue"])
d2g_mid = np.linspace(0.001, 1.3, 100000)
mask = np.where(d2g_mid >= 1, 1, 0)
switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / 0.03))
switch1 = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / 0.01))
ax.plot(d2g_mid, mask, ls='-', color=colors[0], label=r'$n \to 0$')
ax.plot(d2g_mid, switch1, color=colors[1], ls='--', label=r'$n = 0.01$')
ax.plot(d2g_mid, switch, ls='-.', color=colors[2], label=r'$n = 0.03$')

ax.tick_params(axis='both', which='both')
# ax.grid()
ax.legend()
ax.set_xlim(0.3, 1.3)
ax.set_xlim(0.8, 1.2)
#ax.set_ylim(1e-15, 5)
ax.set_xlabel(r"Midplane dust-to-gas ratio ($\rho_{\rm d}/\rho_{\rm g}$)")
ax.set_ylabel(r"Probability ($\mathcal{P}_{\rm pf})$")
filename = localDir + '/figplots/pf_prob_plot'
plt.tight_layout()
plt.savefig(filename+'.png', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.savefig(filename+'.eps', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.show()

