import os
import tempfile
import subprocess
import traceback
import logging
import shutil
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import ticker
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams
M_earth = 5.9722e24 * 1e3  # [g]

width_inches = 6 # inches
golden_ratio = 3/4
height_inches = golden_ratio * width_inches
fontsize = 14


def movieBump(z, dataDir, dir, sd=True, fps=35, still=None):
    """
    Create a movie of the dustpy files residing in the directory `dir`.

    The movie needs `ffmpeg` to be installed and be present in the path.

    You can set other values, e.g. by setting the defaults before calling
    the function, or by using a context:
    #
    # >>> with plt.style.context({'figure.dpi':150}):
    # >>>     movie('output/movies/')

    Arguments:
    ----------

    dir : str
        directory where the data files reside

    Keywords:
    ---------

    files : str or list
        file name (wildcards allowed) or list of file names in `dir` to read

    showinfall : bool
        whether to show infall region

    size_limits : bool
        whether to show the fragmentation and drift barrier

    showst1 : bool
        whether to show the St=1 line

    fps : int
        frame rate passed to ffmpeg

    yaxis : str
        'mass', 'size', or 'stokes'

    xlim : None | list
        if given, use those as xlimits

    ylim : None | list
        if given, use those as ylimits

    zlim : None | list
        if given, use those as z-limits. unlike x- and ylim, they are in
        log-space, so give for example [-9,1] instead of [1e-9,1e1].

    xlog, ylog : bool
        set to False to get a linear x or y axis scale

    still : None | int
        if not `None`, then just show a y of the given snapshot

    Output:
    -------

    movie is created in the current directory
    """
    localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'
    w.datadir = dataDir

    # Overlay seaborn's styling with personal adjustments
    plt.style.use('seaborn-paper')
    plt.style.use('tex')
    plt.rcParams["figure.figsize"] = width_inches, height_inches

    # define name of movie file and create a temporary directory for the images
    if sd:
        moviename = localDir + '/simplots/movies/r' + str(z) + 'detailed'
    else:
        moviename = localDir + '/simplots/movies/d' + str(z)
    if still is None:
        print('creating movie {}.mp4'.format(moviename))
        tempdir = tempfile.mkdtemp(dir=os.curdir)
        i_ini = 0
    else:
        i_ini = still

    # Time data
    t = w.read.sequence('t') / c.year
    Nt = t.shape[0]
    tMyr = t / 1e6
    tMyrEnd = tMyr[Nt - 1]
    print("tMyrEnd = ", tMyrEnd)

    # Mass data
    m = w.read.sequence('grid.m')  # Mass grid field [g]
    Nm = m.shape[1]  # Number of mass bins
    A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]  # Grid constant
    dm = 2. * (A - 1.) / (A + 1.)  # mass bin width

    # Radial data
    r = w.read.sequence('grid.r')  # Radial grid cell centers [cm]
    R = r / c.au  # Radial grid cell centers [AU]
    Nr = R.shape[1]
    rInt = w.read.sequence('grid.ri')  # Radial grid cell interfaces [cm]

    # Misc data
    alpha = w.read.sequence('gas.alpha')
    d2g = w.read.sequence('dust.eps')

    # Planetesimal information
    SigmaPlan = w.read.sequence('planetesimals.Sigma')
    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
    PlanMass = w.read.sequence('planetesimals.M') / M_earth
    PlanDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / c.M_sun
    PlanDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / M_earth
    print("Mass of final planetesimal disc mass in Earth masses: %.10f" % PlanDiskMassEarth[-1])

    # Dust information
    SigmaDust = w.read.sequence('dust.Sigma')
    SigmaDustTot = np.sum(SigmaDust, axis=2)
    DustDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / c.M_sun
    DustDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :],
                               axis=1) / M_earth
    print("Initial dust disc mass (Earths): ", DustDiskMassEarth[0])
    print("Final dust disc mass (Earths): ", DustDiskMassEarth[-1])
    SigmaDustDist = SigmaDust / dm
    particleSize = w.read.sequence('dust.a')  # Particle size field [cm]

    # Gas information
    SigmaGas = w.read.sequence('gas.Sigma')
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
    GasDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / M_earth
    SigmaGasDist = SigmaGas / dm

    # Dust-to-gas ratio
    rho_gas = w.read.sequence('gas.rho')
    rho_dust = w.read.sequence('dust.rho').sum(-1)
    d2g_mid = rho_dust / rho_gas

    if sd is False:
        data = w.read.output(Nt - 1)
        filename = "data"
        extension = "hdf5"
        PlanMassEarth = data.planetesimals.M[-1] / M_earth
        eps = data.dust.eps

        cs = data.gas.cs
        delta = data.dust.delta.turb
        OmegaK = data.grid.OmegaK
        St = data.dust.St
        vK = OmegaK * r
        vFrag = data.dust.v.frag

        # Fix indices if necessary
        im = Nm - 1
        ir = Nr - 1
        it = Nt - 1

        # Transformation of distribution
        a = np.mean(m[..., 1:] / m[..., :-1], axis=-1)
        dm = 2. * (a - 1.) / (a + 1.)
        sigmaD = SigmaDust[..., :] / dm

        # Calculating limits

        # Fragmentation limit
        b = vFrag ** 2 / (delta * cs ** 2)
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings(
                'ignore',
                r'invalid value encountered in sqrt')
            St_fr = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b ** 2))

        # Drift limit
        p = SigmaGas * OmegaK * cs

        _f = interp1d(np.log10(r), np.log10(p), fill_value='extrapolate')
        pi = 10. ** _f(np.log10(rInt))
        gamma = np.abs(r / p * np.diff(pi) / np.diff(rInt))
        St_dr = eps / gamma * (vK / cs) ** 2

        # Get limits
        sd_max = np.ceil(np.log10(sigmaD.max()))
        sg_max = np.ceil(np.log10(SigmaGas.max()))

    # getting the data
    [alpha, amplitude, velocity] = getJobParams(z)
    fig = plt.figure()
    gs = gridspec.GridSpec(nrows=1, ncols=1) #, width_ratios=[20,1])

    ptot = f"{PlanDiskMassEarth[Nt - 1]:.1f}"
    textstr = "Plan Mass: " + str(ptot) + " Earths"
    stationary = False
    if z >= 150:
        stationary = True
    if stationary:
        basetitle = r"$\alpha$" + "={a}, A={A}, r$_p$={v}au, ".format(a=alpha, A=amplitude, v=velocity)
    else:
        basetitle = r"$\alpha$" + "={a}, A={A}, v={v}%, ".format(a=alpha, A=amplitude, v=velocity)
    if still is None:
        # the update/image creation loop
        if z <= 201:
            [alpha0, amplitude, position] = getJobParams(z)
        else:
            position = 90
        iguess = np.argmin(abs(R[-1] - position))
        for it in range(Nt):
            igap = np.argmin(SigmaGas[it, 0:iguess + 10])
            iguess = igap
            ipeak = np.argmax(SigmaDustTot[it, igap:igap + 35]) + igap
            dist = ipeak - igap
            istart = int(ipeak - 0.5 * dist)
            iend = int(ipeak + 0.5 * dist)

            ax0 = fig.add_subplot(gs[0, 0])
            titlestr = basetitle + '{t:.2f} Myr'.format(t=tMyr[it])
            ax0.set_title(titlestr, fontdict={'fontsize': fontsize})
            ax0.text(0.05, 0.92, textstr, transform=ax0.transAxes, fontsize=10)
            ax0.set_xlabel("Distance from star [au]")
            fig.tight_layout()

            if sd is True:
                ymin = 1e-6
                ymax = 1e4
                ax0.set_ylim(ymin, ymax)
                ax0.set_ylabel("Surface Density [g/cm²]")
                ax0.loglog(R[it, ...], SigmaDustTot[it, ...], label="Dust")
                ax0.loglog(R[it, ...], SigmaGas[it, ...], label="Gas")
                ax0.loglog(R[it, ...], SigmaPlan[it, ...], label="Planetesimals")
                ax0.loglog(R[it, ...], d2g[it, ...], label="Dust-to-gas ratio")
                ax0.vlines(R[it, ipeak], ymin, ymax, 'r')
                ax0.vlines(R[it, istart], ymin, ymax, 'gray', '--')
                ax0.vlines(R[it, iend], ymin, ymax, 'gray', '--')
                # ax0.legend(loc='upper right')
            else:
                pltcmap = ax0.contourf(r / c.au,
                                          m,
                                          np.log10(sigmaD[it]),
                                          levels=np.linspace(sd_max - 6, sd_max, 7),
                                          cmap="magma",
                                          extend="both"
                                          )
                ax0.contour(r / c.au,
                               m,
                               St[it],
                               levels=[1.],
                               colors="white",
                               linewidths=2
                               )
                ax0.contour(r / c.au,
                               m,
                               (St - St_dr[..., None])[it],
                               levels=[0.],
                               colors="C2",
                               linewidths=1
                               )
                ax0.contour(r / c.au,
                               m,
                               (St - St_fr[..., None])[it],
                               levels=[0.],
                               colors="C0",
                               linewidths=1
                               )

                ax0.axhline(m[im], color="#AAAAAA", lw=1, ls="--")
                ax0.axvline(r[ir] / c.au, color="#AAAAAA", lw=1, ls="--")
                cbarcmap = plt.colorbar(pltcmap, ax=ax0)
                cbarcmap.ax.set_ylabel("$\log\ \sigma$ [g/cm²]")
                ax0.set_xlim(r[0] / c.au, r[-1] / c.au)
                ax0.set_ylim(m[0], m[-1])
                ax0.set_xscale("log")
                ax0.set_yscale("log")
                ax0.set_xlabel("Distance from star [au]")
                ax0.set_ylabel("Particle mass [g]")

            # save image
            fig.savefig(os.path.join(tempdir,'img_{:03d}.png'.format(it))) #bbox_inches="tight") #, pad_inches=0.05)
            plt.cla()
            plt.clf()

        # create movie
        ret = subprocess.call(['ffmpeg', '-i', os.path.join(tempdir, 'img_%03d.png'), '-c:v', 'libx264',
                               '-crf', '15', '-maxrate', '400k', '-pix_fmt', 'yuv420p', '-r', str(fps),
                               '-bufsize', '1835k', '-t', '00:00:05', moviename + '.mp4'])
        shutil.rmtree(tempdir, ignore_errors=True)
        # delete_imgs = ''
        # try:
        #     ret = subprocess.call(['ffmpeg', '-i', os.path.join(tempdir,'img_%03d.png'), '-c:v', 'libx264',
        #                            '-crf', '15', '-maxrate', '400k', '-pix_fmt', 'yuv420p', '-r', str(fps),
        #                            '-bufsize', '1835k', '-t', '00:00:05', moviename + '.mp4'])
        #     if ret == 0:
        #         delete_imgs = 'y'
        #         print('Movie created successfully')
        # except Exception as e:
        #
        #     # display the error
        #
        #     logging.error(traceback.format_exc())
        #     delete_imgs = ''
        #
        # finally:
        #
        #     # ask if images should be kept
        #
        #     while delete_imgs.lower() not in ['y', 'n']:
        #         delete_imgs = input("Error during movie creation, delete images [y/n]").lower()
        #
        #     if delete_imgs == 'y':
        #         print('Deleting images.')
        #         shutil.rmtree(tempdir, ignore_errors=True)
        #     else:
        #         print('Keeping images. They can be found in folder \'{}\''.format(tempdir))


def num2tex(n, x=2, y=2):
    r"""This function turns a real number into a tex-string numbers >10^x and <10^-x
    are returned as e.g., $1.23 \times 10^{5}$, otherwise as e.g., $1234$.
    Unnecessary digit in the tex-string are removed

    Arguments:
    ----------
    n = number to be converted to tex string

    Keywords:
    ---------
    x : int
    :    threshold exponent

    y : int
    :    number of digits

    Example:
    --------
    >>> num2tex([3e3,3e4,3e5],5,1)
    '$3000.0$, $30000.0$, $3\times 10^{5}$'

    """
    from numpy import array, log10
    s = None
    for i in array(n, ndmin=1):
        if i == 0:
            t = r'$0$'
        else:
            if log10(i) > x or log10(i) < -x:
                t = ('%2.' + str(y) + 'e') % i
                t = t[0:t.find('e')]
                t = r'$' + t + r' \times 10^{%i}$' % round(log10(i / float(t)))
            else:
                t = ('%' + str(y) + '.' + str(y) + 'f') % i
                t = r'$' + t + '$'
        #
        # some corrections
        #
        if y == 0:
            nz = ''
        else:
            nz = '.' + str(0).zfill(y)
        t = t.replace('1' + nz + ' \times ', '')
        t = t.replace(nz + ' ', '')
        #
        # we don't need 1\times 10^{x}, just 10^{x}
        #
        t = t.replace(r'$1\times', '$')
        #
        # if there are more input numbers, attache them with a comma
        #
        if s is not None:
            s = s + ', ' + t
        else:
            s = t
    return s