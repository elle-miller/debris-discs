import os
import tempfile
import subprocess
import traceback
import logging
import shutil

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import ticker
from dustpy import hdf5writer as w
from dustpy import readdump
from dustpy import constants as c
from jobInfo import getJobParams


def movieBump(z,dir, files="data*.hdf5", showinfall=False, size_limits=True, fps=20,
          yaxis='mass', ylim=None, showst1=False, xlim=None, zlim=None,
          xlog=True, ylog=True,
          still=None):
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
    w.datadir = localDir + '/sims/' + str(z)

    # define name of movie file and create a temporary directory for the images

    moviename = os.path.normpath(dir).split(os.path.sep)[-1]+"sd"
    if still is None:
        print('creating movie {}.mp4'.format(moviename))
        tempdir = tempfile.mkdtemp(dir=os.curdir)
        i_ini = 0
    else:
        i_ini = still

    # getting the data
    [alpha, amplitude, velocity] = getJobParams(z)
    _r = getSequence("grid/r", files)
    r = _r / c.AU
    t = getSequence("t", files) / c.yr
    tMyr = t / 1e6
    Nt = t.shape[0]
    tMyrEnd = tMyr[Nt - 1]
    d2g = getSequence("dust/dust2gasRatio", files)
    rInt = getSequence("grid/rInt", files)
    # Read the mass grid m,  has Nt rows and Nm columns
    m = getSequence("grid/m", files)
    Nm = m.shape[1]
    M_earth = 5.9722e24 * 1e3  # grams

    # Compute the grid constant and the mass bin width
    A = np.mean(m[:, 1:] / m[:, :-1], axis=1)[..., None, None]
    dm = 2. * (A - 1.) / (A + 1.)
    # The DustPy units are such that the mass is integrated over a mass bin. That
    # means numerically we can sum up the mass dimension to get to dust surface density
    # SigmaDust and SigmaDustTot has Nt rows and Nr cols
    # Divide by dm to get the logarithmic density distribution
    SigmaDust = getSequence("dust/Sigma", files)
    SigmaDustTot = np.sum(SigmaDust, axis=2)
    DustDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaDustTot[:, :], axis=1) / c.M_sun
    SigmaDustDist = SigmaDust / dm  # 4 x 100
    particleSize = getSequence("dust/a", files)  # Nt x Nr = 4 x 100

    # Gas information, SigmaGas_m has Nt rows and Nr columns
    SigmaGas = getSequence("gas/Sigma", files)
    SigmaGasTot = np.sum(SigmaGas, axis=-1)
    GasDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaGas[:, :], axis=1) / c.M_sun
    SigmaGasDist = SigmaGas / dm

    # Planetesimal information
    SigmaPlan = getSequence("dust/SigmaPlan", files)
    SigmaPlanTot = np.sum(SigmaPlan, axis=-1)
    PlanDiskMass = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / c.M_sun
    PlanDiskMassEarth = np.sum(np.pi * (rInt[:, 1:] ** 2. - rInt[:, :-1] ** 2.) * SigmaPlan[:, :], axis=1) / M_earth


    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1])
    if still is None:
        # the update/image creation loop

        for it in range(len(files)):

            ax0 = fig.add_subplot(gs[0, 0])
            ax0.set_ylim(1.e-6, 1.e4)
            ax0.set_xlabel("Distance from star [AU]")
            ax0.set_ylabel("Surface Density [g/cm²]")
            ti = ax0.set_title('time')
            fig.tight_layout()

            ax0.loglog(r[-1, ...], SigmaDustTot[it, ...], label="Dust")
            ax0.loglog(r[-1, ...], SigmaGas[it, ...], label="Gas")
            ax0.loglog(r[-1, ...], SigmaPlan[it, ...], label="Planetesimals")
            ax0.loglog(r[-1, ...], d2g[it, ...], label="Dust-to-gas ratio")
            ax0.legend()

            # update title text

            ti.set_text(r"$\alpha$" + "={a}, A={A}, v={v}%, ".format(a=alpha, A=amplitude, v=velocity) + 't = {} yr'.format(num2tex(t[it])))

            # save image

            fig.savefig(os.path.join(tempdir,'img_{:03d}.png'.format(it)))
            plt.cla()
            plt.clf()

        # create movie
        delete_imgs = ''
        try:
            ret = subprocess.call(['ffmpeg', '-i', os.path.join(tempdir,'img_%03d.png'), '-c:v', 'libx264',
                                   '-crf', '15', '-maxrate', '400k', '-pix_fmt', 'yuv420p', '-r', str(fps),
                                   '-bufsize', '1835k', moviename + '.mp4'])
            if ret == 0:
                delete_imgs = 'y'
                print('Movie created successfully')
        except Exception as e:

            # display the error

            logging.error(traceback.format_exc())
            delete_imgs = ''

        finally:

            # ask if images should be kept

            while delete_imgs.lower() not in ['y', 'n']:
                delete_imgs = input("Error during movie creation, delete images [y/n]").lower()

            if delete_imgs == 'y':
                print('Deleting images.')
                shutil.rmtree(tempdir, ignore_errors=True)
            else:
                print('Keeping images. They can be found in folder \'{}\''.format(tempdir))


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