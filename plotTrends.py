import dustpy
from dustpy import constants as c
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter, StrMethodFormatter, NullFormatter
import matplotlib.ticker as ticker
from dustpy import hdf5writer as w
from plottingFunctions import *

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = getcwd()
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'
outputDir = localDir + '/figplots/'

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14

plt.style.use('seaborn-paper')
plt.style.use('tex')
plt.rcParams["figure.figsize"] = width_inches, height_inches

# # Figure stuff
markersize = 12
fontsize = 14
labelsize = 14
length = 3
width = 0.5

# Center, ring width and fractional width data for moving bump
floor = 0
velocity = [10, 30, 100]
centerList310 = [108.1, 92, 51.2]
widthList310 = [12, 31.9, 64.3]
fracList310 = [0.11, 0.35, 1.26]
widthList410 = [4.2, 6.3, 13.1]
centerList410 = [114.1, 113.1, 109.2]
fracList410 = [0.04, 0.06, 0.13]

# Ring start and end information for moving bump
vel2 = [10, 30, 50, 75, 100, 300]
vel3 = [10, 30, 100, 300]
vel = [10, 30, 100, 300]
start310 = [114.12, 98.464, 52.59, 19.425]
end310 = [98.464, 54.57, 14.46, 15.568]
start330 = [122.86, 102.16, 84.956, 65.622, 54.57, 16.76]
end330 = [102.16, 56.62, 15, 15, 14.46, 16.1535]
start410 = [122.857, 118.4077, 118.4077, 118.4077]
end410 = [114.12, 109.986, 98.4639, 54.569]
start430 = [127.474, 127.474, 127.474, 127.474]
end430 = [122.857, 118.4077, 102.164, 56.62]

# plotPositions3 = True
# if plotPositions3:
#     vel4 = [0, 10, 30, 100, 300]
#     start310 = [114.12, 98.464, 52.59, 19.425]
#     end310 = [98.464, 54.57, 14.46, 15.568]
#     start330 = [122.86, 102.16, 54.57, 16.76]
#     end330 = [102.16, 56.62, 14.46, 16.1535]
#     start410 = [118.4077, 122.857, 118.4077, 118.4077, 118.4077]
#     end410 = [114.12, 114.12, 109.986, 98.4639, 54.569]
#     start430 = [127.474, 127.474, 127.474, 127.474, 127.474]
#     end430 = [122.857, 122.857, 118.4077, 102.164, 56.62]
#     fig, ax = plt.subplots()
#     ax.plot(vel, start330, markersize=markersize, marker='*', ls='--', color="C0", label=r"$\alpha_0$ = 1e-3")
#     ax.plot(vel, end330, markersize=markersize, marker='*', ls='--', color="C0")
#     ax.plot(vel4, start430, markersize=markersize, marker='.', ls='-.', color="C3", label=r"$\alpha_0$ = 1e-4")
#     ax.plot(vel4, end430, markersize=markersize, marker='.', ls='-.', color="C3")
#     ax.tick_params(axis='y', which='both', length=length, width=width, labelsize=labelsize)
#     ax.tick_params(axis='x', which='both', length=length, width=width, labelsize=labelsize)
#     ax.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)
#     ax.set_ylabel("Final ring location [AU]", fontsize=fontsize)
#     ax.legend(fontsize=fontsize - 4, loc='lower left')
#     ax.set_ylim(0, 140)
#     ax.set_title("A = 30 bump at 90 au evolved for 10 Myr", fontsize=fontsize)
#     fig.tight_layout()
#     plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/ring10.png', format='png',
#                 dpi=600)
#     plt.show()

plot34 = True
if plot34:

    fig, ax1 = plt.subplots()

    # Data
    center = 3
    n = 3

    if center == 1:
        ax1.set_ylabel("Ring width [au]", fontsize=fontsize)
        loc = 'upper left'
        figname = 'width34'
        lns3 = ax1.semilogx(velocity[0:n], widthList310[0:n], markersize=markersize, marker='^', ls='--', color="C2",
                            label=r"$\alpha$ = 1e-3")
        lns4 = ax1.semilogx(velocity[0:n], widthList410[0:n], markersize=markersize, marker='*', ls='--', color="C3",
                            label=r"$\alpha$ = 1e-4")
        lns = lns3 + lns4
    if center == 2:
        ax1.set_ylabel("Ring center [au]", fontsize=fontsize)

        loc = 'lower left'

        figname = 'center34'
        lns3 = ax1.plot(velocity[0:n], centerList310[0:n], markersize=markersize, marker='o', ls='--', color="C0",
                        label=r"$\alpha$ = 1e-3")
        lns4 = ax1.plot(velocity[0:n], centerList410[0:n], markersize=markersize, marker='v', ls='-.', color="C1",
                        label=r"$\alpha$ = 1e-4")
        lns = lns3 + lns4
    if center == 3:
        ax1.set_ylabel("Fractional width", fontsize=fontsize)
        loc = 'upper left'

        figname = 'frac34'
        lns3 = ax1.plot(velocity[0:n], fracList310[0:n], markersize=markersize, marker='p', ls='--', color="C6",
                        label=r"$\alpha$ = 1e-3")
        lns4 = ax1.plot(velocity[0:n], fracList410[0:n], markersize=markersize, marker='>', ls='-.', color="C7",
                        label=r"$\alpha$ = 1e-4")
        lns = lns3 + lns4

    # Legend

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=loc, fontsize=fontsize)

    # Axis
    titlestr = "Ring from moving A=10 bump evolved for 10 Myr"
    ax1.set_title(titlestr, fontsize=fontsize)
    ax1.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)

    # Tickers
    ax1.tick_params(axis='y', which='both', length=length, width=width, labelsize=labelsize)
    ax1.tick_params(axis='x', which='both', length=length, width=width, labelsize=labelsize)

    # Formatting
    fig.tight_layout()
    plt.savefig(outputDir + figname + '.png', format='png',
                dpi=600)
    plt.show()

plotRingChars4 = True
if plotRingChars4:

    fig, ax1 = plt.subplots()

    # Data
    center = 2
    alpha = 3
    if alpha == 3:
        n = 3
    else:
        n = 4

    if center == 1:
        ax1.set_ylabel("Ring Width [AU]", fontsize=fontsize)
        loc = 'upper left'
        figname = 'widthVelocity4'
        lns3 = ax1.semilogx(velocity[0:n], widthList410[0:n], markersize=markersize, marker='^', ls='None', color="C2",
                            label="A = 10")
        lns4 = ax1.semilogx(velocity[0:n], widthList430[0:n], markersize=markersize, marker='*', ls='None', color="C3",
                            label="A = 30")
        lns = lns3 + lns4
    if center == 2:
        ax1.set_ylabel("Ring Center [AU]", fontsize=fontsize)
        if alpha == 3:
            loc = 'upper right'

            figname = 'centerVelocity3'
            lns3 = ax1.plot(velocity[0:n], centerList310[0:n], markersize=markersize, marker='o', ls='--', color="C0",
                            label="A = 10")
            lns4 = ax1.plot(velocity[0:n], centerList330[0:n], markersize=markersize, marker='v', ls='-.', color="C1",
                            label="A = 30")
        else:
            figname = 'centerVelocity4'
            loc = 'lower left'
            lns3 = ax1.plot(velocity[0:n], centerList410[0:n], markersize=markersize, marker='o', ls='--', color="C0",
                            label="A = 10")
            lns4 = ax1.plot(velocity[0:n], centerList430[0:n], markersize=markersize, marker='v', ls='-.', color="C1",
                            label="A = 30")
        lns = lns3 + lns4
    if center == 3:
        ax1.set_ylabel("Fractional Width", fontsize=fontsize)
        loc = 'upper left'
        if alpha == 4:
            figname = 'fracVelocity4'
            lns3 = ax1.plot(velocity[0:n], fracList410[0:n], markersize=markersize, marker='p', ls='--', color="C6",
                            label="A = 10")
            lns4 = ax1.plot(velocity[0:n], fracList430[0:n], markersize=markersize, marker='>', ls='-.', color="C7",
                            label="A = 30")
        else:
            figname = 'fracVelocity3'
            lns3 = ax1.plot(velocity[0:n], fracList310[0:n], markersize=markersize, marker='p', ls='--', color="C6",
                            label="A = 10")
            lns4 = ax1.plot(velocity[0:n], fracList330[0:n], markersize=markersize, marker='>', ls='-.', color="C7",
                            label="A = 30")
        lns = lns3 + lns4

    # Legend

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=loc, fontsize=fontsize - 2)

    # Axis
    titlestr = "Ring from " + r"$\alpha_0$ = 1e-" + str(alpha) + " bump at 90 au evolved for 10 Myr"
    ax1.set_title(titlestr, fontsize=fontsize)
    ax1.set_xlabel("Bump velocity as \% of nominal", fontsize=fontsize)

    # Tickers
    ax1.tick_params(axis='y', which='both', length=length, width=width, labelsize=labelsize)
    ax1.tick_params(axis='x', which='both', length=length, width=width, labelsize=labelsize)

    # Formatting
    fig.tight_layout()
    plt.savefig(outputDir + figname + '.png', format='png',
                dpi=600)
    plt.show()

# if plotVels2:
#
# 	fig, ax1 = plt.subplots()
# 	ax2 = ax1.twinx()
#
# 	# Data
# 	n=9
# 	lns1 = ax1.loglog(velocity[-n:],fracList[-n:],ls='-',color="C2",label="Fractional Width")
# 	lns2 = ax2.loglog(velocity[-n:], widthList[-n:],ls='--',color="C0",label="Ring Width")
# 	lns3 = ax2.loglog(velocity[-n:],centerList[-n:],ls='-.',color="C1",label="Center Position")
#
# 	# Legend
# 	lns = lns1+lns2+lns3
# 	labs = [l.get_label() for l in lns]
# 	ax1.legend(lns, labs, loc='lower right',fontsize=fontsize)
#
# 	# Axis
# 	ax1.set_xlim(min(velocity[-n:]),max(velocity[-n:]))
# 	ax1.set_xlabel("Planetary Velocity [AU/Myr]",fontsize=fontsize)
# 	ax1.set_ylabel("Fractional Width",fontsize=fontsize)
# 	ax2.set_ylabel("Distance [AU]",fontsize=fontsize)
#
# 	# Tickers
# 	ax1.tick_params(axis='y',which='both',length=length,width=width,labelsize=labelsize)
# 	ax1.tick_params(axis='x',which='both',length=length,width=width,labelsize=labelsize)
# 	ax2.tick_params(axis='y',which='both',length=length,width=width,labelsize=labelsize)
# 	formatter = ticker.ScalarFormatter()
# 	formatter.set_scientific(False)
# 	for axis in [ax1.xaxis, ax1.yaxis,ax2.yaxis]:
# 		axis.set_minor_formatter(ScalarFormatter())
# 		axis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
# 	ax1.yaxis.set_minor_formatter(StrMethodFormatter('{x:.2f}'))
# 	for label in ax2.get_yaxis().get_ticklabels(which='both')[::2]:
#     		label.set_visible(False)
# 	for label in ax1.get_yaxis().get_ticklabels(which='both')[::2]:
# 		label.set_visible(False)
#
# 	# Formatting
# 	fig.tight_layout()
# 	plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/velocity_rings_4.eps', format='eps',dpi=1200)
# 	plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/velocity_rings_4.png', format='png',dpi=1200)
# 	plt.show()
#
# if plotVels1:
#
# 	fig, ax1 = plt.subplots()
# 	ax2 = ax1.twinx()
#
# 	# Data
# 	lns1 = ax1.plot(velocity2,fracList2,ls='-',color="C2",label="Fractional Width")
# 	lns2 = ax2.plot(velocity2, widthList2,ls='--',color="C0",label="Ring Width")
# 	lns3 = ax2.plot(velocity2[0:5], centerList2[0:5],ls='-.',color="C1",label="Center Position")
#
# 	# Legend
# 	lns = lns1+lns2+lns3
# 	labs = [l.get_label() for l in lns]
# 	ax1.legend(lns, labs, loc='upper right',fontsize=fontsize)
#
# 	# Axis
# 	ax1.set_xlim(min(velocity2),max(velocity2))
# 	#ax2.set_ylim(0,max(centerList2))
# 	ax1.set_xlabel("Planetary Velocity [AU/Myr]",fontsize=fontsize)
# 	ax1.set_ylabel("Fractional Width",fontsize=fontsize)
# 	ax2.set_ylabel("Distance [AU]",fontsize=fontsize)
#
# 	# Tickers
# 	ax1.tick_params(axis='y',which='both',length=length,width=width,labelsize=labelsize)
# 	ax1.tick_params(axis='x',which='both',length=length,width=width,labelsize=labelsize)
# 	ax2.tick_params(axis='y',which='both',length=length,width=width,labelsize=labelsize)
# 	ax1.yaxis.set_minor_formatter(ScalarFormatter())
# 	ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
# 	formatter = ticker.ScalarFormatter()
# 	formatter.set_scientific(False)
# 	for axis in [ax1.xaxis, ax2.yaxis]:
# 		axis.set_minor_formatter(ScalarFormatter())
# 		axis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
# 	for label in ax2.get_yaxis().get_ticklabels(which='both')[::2]:
#     		label.set_visible(False)
# 	for label in ax1.get_yaxis().get_ticklabels(which='both')[::2]:
# 		label.set_visible(False)
#
# 	# Formatting
# 	fig.tight_layout()
# 	plt.savefig('/media/elle/Seagate Backup Plus Drive/2020/mpia/mpia/paperplots/velocity_rings_3.png', format='png',dpi=1200)
# 	plt.show()