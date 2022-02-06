from dustpy import plot
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os, os.path
from plottingFunctions import *

# cd /mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs/


# this can be run from where the data is stored. You can do $ python getplanmasses.py -z 100 -a 1 -t 1 to print
# out the plan masses over time or just the final one in a text file

# NOTE: I think this file is redundant as DustPy v2 has a way of reading Mass directly.

# Global settings
M_earth = 5.9722e24 * 1e3  # [g]

localDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'

outputDir = localDir + '/plandata/'


def main(args):
    z = args.z
    w.datadir = getDataDir(z)
    t = w.read.sequence('t') / (c.year * 1e6)
    Nt = t.shape[0]
    num = len(w.read.listfiles())
    [alpha, amplitude, position] = getJobParams(z)

    planmass = w.read.sequence('planetesimals.M') / M_earth
    firstcutoff = 0.9*planmass[-1]
    cutoffmass = 0.95*planmass[-1]
    print("Final plan mass = ", planmass[-1], " Earths")
    first = False
    for it in range(Nt):
        if planmass[it] >= firstcutoff and not first:
            print("Time at 0.90M = ", t[it], " Myrs")
            first = True
        if planmass[it] >= cutoffmass:
            print("Time at 0.95M  = ", t[it], " Myrs")
            break

    # If we want to add it to masterPlanInfo.txt, only need the final output
    if args.writeToAll:

        filename = 'masterPlanInfo.txt'
        print("writing data to " + filename)
        try:
            data = w.read.output(num - 1)
        except:
            data = w.read.output(num - 2)

        t = (data.t / c.year * 1e-6)
        #t = w.read.sequence('t') / c.year * 1e-6
        #pm = w.read.sequence('planetesimals.M') / M_earth
        pm = data.planetesimals.M / M_earth
        text = str(z) + "  " + "{a:6}  {A:2}  {p:3}  {t:4.2f}  {pm:6.2f}".format(a=alpha, A=amplitude, p=position,
                                                                                 t=t[0], pm=pm[0])
        print(text)
        append_new_line(filename, text)

    # Write the evolution of one particular file
    if args.writeTime:
        filename = outputDir + str(z) + 'planinfo.txt'
        pm = w.read.sequence('planetesimals.M') / M_earth
        t = w.read.sequence('t') / c.year * 1e-6

        append_new_line(filename, 'z    alpha  A  pos[au]')
        append_new_line(filename, str(z) + "  " + "{a:6}  {A:2}  {p:3}".format(a=alpha, A=amplitude, p=position))
        append_new_line(filename, 't[Myr] planMass[Earths]')
        for n in range(num):
            text = "{t:4.2f}  {pm:15.10f}".format(t=t[n], pm=pm[n])
            append_new_line(filename, text)


#    if ringWidth1[Nt - 1] == 0:
#     textstr = ""
#    elif ringWidth2[Nt - 1] == 0:
#     textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1
#    else:
#     textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1 + "\n" + "2: c=" + c2 + " w=" + w2 + " AU, f=" + f2

#    # Create strings for plots
#    [alpha, amplitude, position] = getJobParams(z)
#    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
#    textstr = textstr + "\nPlan Disc Mass: " + str(ptot) + " Earths"
#    titlestr = str(z) + ": " + r"$\alpha$" + "={a}, A={A}, $r_p$={p}AU @ {t:.2f} Myr".format(a=alpha, A=amplitude,
#                                                                                                  p=position,
#                                                                                                  t=tMyrEnd)
#    	
#    	
#   
#    	


def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)


#    

#    # Get basic data from files
#    data = hdf5writer.read.all()
#    t = data.t / c.year
#    Nt = t.shape[0]
#    tMyr = t / 1e6
#    tMyrEnd = tMyr[Nt - 1]
#    print("tMyrEnd = ", tMyrEnd)
#    d2g = data.dust.eps
#    rInt = data.grid.ri  # Radial grid cell interfaces [cm]
#    m = data.grid.m  # Mass grid field [g]
#    Nm = m.shape[1]  # Number of mass bins
#    A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]  # Grid constant
#    dm = 2. * (A - 1.) / (A + 1.)  # mass bin width
#    r = data.grid.r  # Radial grid cell centers [cm]
#    R = r / c.au  # Radial grid cell centers [AU]
#    Nr = R.shape[1]

#    # Dust information
#    SigmaDust = data.dust.Sigma
#    SigmaDustTot = np.sum(SigmaDust, axis=2)
#    DustDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / c.M_sun
#    DustDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / M_earth
#    print("Initial dust disc mass (Earths): ", DustDiskMassEarth[0])
#    print("Final dust disc mass (Earths): ", DustDiskMassEarth[-1])

#    SigmaDustDist = SigmaDust / dm
#    particleSize = data.dust.a  # Particle size field [cm]

#    # Gas information
#    SigmaGas = data.gas.Sigma
#    SigmaGasTot = np.sum(SigmaGas, axis=-1)
#    GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
#    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
#    SigmaGasDist = SigmaGas / dm

#    # Planetesimal information
#    SigmaPlan = data.planetesimals.Sigma
#    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
#    PlanMass = data.planetesimals.M / M_earth
#    PlanDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / c.M_sun
#    PlanDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / M_earth
#    print("Mass of final planetesimal disc mass in Earth masses: %.10f" % PlanDiskMassEarth[-1])
#    
#    
##     first instance of planetesimals
#     # for j in range(Nt):
#     #    for i in range(Nr):
#     #        if SigmaPlan[j,i]>1e-99:
#     #            print("bingo")
#     #            print("time = %f" % (t[i]/1e6))
#     #            print("pos = %f" % R[j,i])
#     #            print("d2g = %f" % d2g[j,i])
#     # Initialize
#    ringWidth1 = np.zeros(Nt)
#    startRing1 = np.zeros(Nt)
#    endRing1 = np.zeros(Nt)
#    centerRing1 = np.zeros(Nt)
#    fracWidth1 = np.zeros(Nt)
#    ringWidth2 = np.zeros(Nt)
#    startRing2 = np.zeros(Nt)
#    endRing2 = np.zeros(Nt)
#    centerRing2 = np.zeros(Nt)
#    fracWidth2 = np.zeros(Nt)

#    # Duss in Ring
#    RingDustTot = SigmaDustTot.copy()
#    Ring1DustTot = SigmaDustTot.copy()
#    Ring2DustTot = SigmaDustTot.copy()

#    # For every epoch
#    twoRingsFlag = False  # Bad coding practice, basically once two rings are flagged this forces always two rings
#    formed = False
#    for j in range(Nt):

#     # Reset index positions
#     iStartRing1P = 0
#     iEndRing1P = 0
#     iStartRing2P = 0
#     iEndRing2P = 0
#     index = 0
#     floorVal = 1.0e-90
#     minVal = min(SigmaPlan[j])  # 1e-100


#     # Loop through each radial bin, locating index positions of start and end ring
#     numRings = 0
#     beginRing2 = False
#     for i in SigmaPlan[j]:
#         if (i > minVal) & (iStartRing1P == 0):
#             if formed is False:
#                 formationTimeIndex = j
#                 formed = True
#             iStartRing1P = index
#             numRings = 1
#         elif (beginRing2 == False) & (numRings == 1) & (i < floorVal) & (iStartRing1P != 0):
#             iEndRing1P = index
#             beginRing2 = True
#         elif beginRing2 & (i > minVal) & (iStartRing2P == 0):
#             iStartRing2P = index
#             numRings = 2
#         elif beginRing2 & (i < floorVal) & (iStartRing2P != 0):
#             iEndRing2P = index
#             break
#         index += 1

#     startRing1[j] = 1e-10
#     endRing1[j] = 1e-10
#     startRing2[j] = 1e-10
#     endRing2[j] = 1e-10

#     # Convert these indices to actual values
#     if iStartRing1P != 0:
#         startRing1[j] = rInt[j, iStartRing1P] / c.AU
#         endRing1[j] = rInt[j, iEndRing1P] / c.AU

#     if iStartRing2P != 0:
#         startRing2[j] = rInt[j, iStartRing2P] / c.AU
#         endRing2[j] = rInt[j, iEndRing2P] / c.AU

#     centerRing1[j] = (endRing1[j] + startRing1[j]) / 2
#     ringWidth1[j] = endRing1[j] - startRing1[j]
#     fracWidth1[j] = ringWidth1[j] / centerRing1[j]
#     centerRing2[j] = (endRing2[j] + startRing2[j]) / 2
#     ringWidth2[j] = endRing2[j] - startRing2[j]
#     fracWidth2[j] = ringWidth2[j] / centerRing2[j]

#     # Now set values not in dust ring to 0
#     # fw=fwhm(R[j],SigmaPlan[j])
#     # iStartRingP = (np.abs(R - fw[1])).argmin()
#     # iEndRingP =  (np.abs(R - fw[0])).argmin()

#     # These are the index values of dust calculated using different way, FWHM
#     # fw=fwhm(R[j],SigmaDustTot[j])
#     # iStartRingD = (np.abs(R - fw[1])).argmin()
#     # iEndRingD =  (np.abs(R - fw[0])).argmin()

#     # For this epoch, loop through all radial bins, and turn OFF dust outside ring indices
#     # print(numRings)
#     # print(iStartRing1P)
#     # print(iEndRing1P)
#     # print(iStartRing2P)
#     # print(iEndRing2P)
#     for k in range(Nr):
#         # If there are no rings, turn off everything
#         if (numRings == 0):
#             RingDustTot[j, :] = 0
#             Ring1DustTot[j, :] = 0
#             Ring2DustTot[j, :] = 0
#             break
#         # If there is only one ring and we aint in it
#         elif (numRings == 1) & ((k not in range(iStartRing1P, iEndRing1P + 1))):
#             RingDustTot[j, k] = 0
#             Ring1DustTot[j, k] = 0
#             Ring2DustTot[j, :] = 0
#         # If there are two rings, and we are in first, middle or end section then set to zero
#         elif (numRings == 2) & ((k < iStartRing1P) | (k in range(iEndRing1P + 1, iStartRing2P)) | (k > iEndRing2P)):
#             RingDustTot[j, k] = 0
#             Ring1DustTot[j, k] = 0
#             Ring2DustTot[j, k] = 0
#         # If in first ring, set ring2 to zero
#         elif (numRings == 2) & (k in range(iStartRing1P, iEndRing1P + 1)):
#             Ring2DustTot[j, k] = 0
#         elif (numRings == 2) & (k in range(iStartRing2P, iEndRing2P + 1)):
#             Ring1DustTot[j, k] = 0

#    RingDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * RingDustTot[:, :], axis=1) / c.M_sun
#    Ring1DiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * Ring1DustTot[:, :], axis=1) / c.M_sun
#    Ring2DiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * Ring2DustTot[:, :], axis=1) / c.M_sun
#    # print("Time of first formation: %.2f" % tMyr[formationTimeIndex])
#    print("Ring 1 Start: %.1f" % startRing1[-1])
#    print("Ring 1 End: %.1f" % endRing1[-1])
#    print("Ring 1 Center: %.1f AU" % centerRing1[Nt - 1])
#    print("Ring 1 Width: %.1f AU" % ringWidth1[Nt - 1])
#    print("Ring 1 Fractional Width: %.2f" % fracWidth1[Nt - 1])
#    print("******")
#    for i in range(Nt-1):
#     print("Time = %.2f" % tMyr[i])
#     print("Plan Mass = %.2f" % PlanDiskMassEarth[i])

#    if ringWidth2[Nt - 1] > 0:
#     print("Ring 2 Center: %.1f AU" % centerRing2[Nt - 1])
#     print("Ring 2 Width: %.1f AU" % ringWidth2[Nt - 1])
#     print("Ring 2 Fractional Width: %.2f" % fracWidth2[Nt - 1])

#    c1 = f"{centerRing1[Nt - 1]:.1f}"
#    w1 = f"{ringWidth1[Nt - 1]:.1f}"
#    f1 = f"{fracWidth1[Nt - 1]:.2f}"
#    c2 = f"{centerRing2[Nt - 1]:.0f}"
#    w2 = f"{ringWidth2[Nt - 1]:.0f}"
#    f2 = f"{fracWidth2[Nt - 1]:.2f}"
#    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"

#    if ringWidth1[Nt - 1] == 0:
#     textstr = ""
#    elif ringWidth2[Nt - 1] == 0:
#     textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1
#    else:
#     textstr = "1: c=" + c1 + " w=" + w1 + " AU, f=" + f1 + "\n" + "2: c=" + c2 + " w=" + w2 + " AU, f=" + f2

#    # Create strings for plots
#    [alpha, amplitude, position] = getJobParams(z)
#    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
#    textstr = textstr + "\nPlan Disc Mass: " + str(ptot) + " Earths"
#    titlestr = str(z) + ": " + r"$\alpha$" + "={a}, A={A}, $r_p$={p}AU @ {t:.2f} Myr".format(a=alpha, A=amplitude,
#                                                                                                  p=position,
#                                                                                                  t=tMyrEnd)
#    # Plot the surface density of dust and gas vs the distance from the star
#    if args.plotSDR:
#        fig, ax = plt.subplots()
#        it = 0
#        ax.loglog(R[-1, ...], SigmaDustTot[it, ...], label="Dust")
#        ax.loglog(R[-1, ...], SigmaGas[it, ...], label="Gas")
#        ax.loglog(R[-1, ...], SigmaPlan[-1, ...], label="Planetesimals")
#        ax.loglog(R[-1, ...], d2g[it, ...], label="d2g Ratio")
#        ax.set_ylim(1.e-6, 1.e4)
#        ax.set_xlabel("Distance from star [AU]")
#        ax.set_ylabel("Surface Density [g/cm²]")
#        ax.legend()
#        ax.set_title(titlestr)
#        ax.text(0.05, 0.9, textstr, transform=ax.transAxes, fontsize=10)
#        fig.tight_layout()
#        filename = outputDir + 'sdr/r' + str(args.z) + '.png'
#        plt.savefig(filename, format='png', dpi=600)
#        plt.show()
#        
#     # Time evolution of gas and dust disk mass
#    if args.plotMass:
#         fig02, ax02 = plt.subplots()
#         ax02.loglog(t, GasDiskMassEarth, label="Gas", color="C0")
#         ax02.loglog(t, DustDiskMassEarth, label="Dust", color="C4")
#         if numRings == 2:
#             ax02.loglog(t, RingDiskMass * c.M_sun / M_earth, ls='--', label="Total Ring Dust", color="C1")
#             ax02.loglog(t, Ring1DiskMass * c.M_sun / M_earth, ls='-.', label="Ring 1 Dust", color="C3")
#             ax02.loglog(t, Ring2DiskMass * c.M_sun / M_earth, ls=':', label="Ring 2 Dust", color="C5")
#         else:
#             ax02.loglog(t, RingDiskMass * c.M_sun / M_earth, ls='--', label="Ring Dust", color="C1")
#         ax02.loglog(t, PlanDiskMassEarth, label="Planetesimals", color="C2")
#         xlim0 = t[min(1, len(t) - 1)]
#         xlim1 = t[-1]
#         ax02.set_xlim(xlim0, xlim1)
#         # ylim0 = 10. ** np.floor(np.log10(np.min(np.minimum(GasDiskMass, DustDiskMass, PlanDiskMass))))
#         # ylim1 = 10. ** np.ceil(np.log10(np.max(np.maximum(GasDiskMass, DustDiskMass, PlanDiskMass))))
#         ax02.set_ylim(1e0, 3e5)
#         ax02.legend(loc='lower left')
#         ax02.lineTime = ax02.axvline(t[0], color="C7", zorder=-1, lw=1)
#         ax02.set_title(titlestr)
#         ax02.set_xlabel("Time [yr]")
#         ax02.set_ylabel("Mass [M$_\oplus$]")
#         ax02.grid(b=False)
#         filename = outputDir + 'mass/m' + str(args.z) + '.png'
#         plt.savefig(filename, format='png', dpi=600)
#         plt.show()

# def fwhm(x, y, k=10):
#     """
#     Determine full-width-half-maximum of a peaked set of points, x and y.

#     Assumes that there is only one peak present in the datasset.  The function
#     uses a spline interpolation of order k.
#     """
#     half_max = max(y) / 2
#     s = splrep(x, y - half_max, k=3)
#     roots = sproot(s)

#     if len(roots) > 2:
#         # raise MultiplePeaks("The dataset appears to have multiple peaks, and thus the FWHM can't be determined.")
#         return [0, 0]

#     elif len(roots) < 2:
#         # raise NoPeaksFound("No proper peaks were found in the data set; likely "the dataset is flat (e.g. all zeros).")
#         return [0, 0]

#     else:
#         return [roots[1], roots[0]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-a', action="store", dest="writeToAll", type=int, default=0, help="Append info to textfile")
    parser.add_argument('-t', action="store", dest="writeTime", type=int, default=0, help="Create a text file of info")
    arguments = parser.parse_args()
    main(arguments)
