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


d2g_mid = np.linspace(0.001, 1.25, 100000)
mask = np.where(d2g_mid >= 1, 1, 0)
switch = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / 0.03))
switch1 = 0.5 * (1. + np.tanh((np.log10(d2g_mid)) / 0.01))
ax.semilogy(d2g_mid, mask, ls='--', color="black", linewidth=1, label='step')
ax.semilogy(d2g_mid, switch1, color="C0", ls='-.', label='n = 0.01')
ax.semilogy(d2g_mid, switch, color="C3", label='n = 0.03')


ax.legend()
ax.set_xlim(0.8, 1.25)
ax.set_ylim(1e-3, 5)
ax.set_xlabel("Midplane dust-to-gas ratio")
ax.set_ylabel("Probability")
filename = localDir + '/figplots/pf_prob_semilog_narrow'
plt.tight_layout()
plt.savefig(filename+'.png', format='png', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.savefig(filename+'.eps', format='eps', bbox_inches='tight', pad_inches=0.05, dpi=300)
plt.show()
