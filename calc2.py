import numpy as np
from matplotlib import pyplot as plt
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

# Figure stuff
# colWidth = 244  # pt
# markersize = 12
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


d2g_mid = np.linspace(0.75, 1.25, 1000)
mask = np.where(d2g_mid >= 1, 1, 0)
switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / 0.03))
ax.plot(d2g_mid, switch, color="C3", label='$\epsilon_{pf}$')
ax.plot(d2g_mid, mask, 'k--', linewidth=0.85, label='step')
ax.legend()
ax.set_xlabel("Midplane dust-to-gas ratio")
ax.set_ylabel("Efficiency")
filename = localDir + '/figplots/pf_prob'
plt.tight_layout()
plt.savefig(filename+'.png', format='png', bbox_inches='tight', pad_inches=0, dpi=600)
plt.savefig(filename+'.eps', format='eps', bbox_inches='tight', pad_inches=0, dpi=600)
plt.show()
