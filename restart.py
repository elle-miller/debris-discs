from dustpy import readdump
from os import path
from dustpy import constants as c
import argparse
import numpy as np

slurmDir = '/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs'
localDir = '/media/elle/Seagate Backup Plus Drive/2020/mpia/debris-discs'
localDirNew = '/media/elle/Seagate Expansion Drive/MPIAResults'
slurmDir = '/mnt/beegfs/bachelor/groups/henning/users/miller/debris-discs'


def main(args):
    z = args.z
    if path.exists(localDir + '/sims/' + str(z)):
        dataDir = localDir + '/sims/' + str(z)
    elif path.exists(localDirNew + '/' + str(z)):
        dataDir = localDirNew + '/' + str(z)
    elif path.exists(slurmDir + '/sims/' + str(z)):
        dataDir = slurmDir + '/sims/' + str(z)
    else:
        print("Output directory not found")
        return

    s = readdump(dataDir + '/frame.dmp')
    s.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', action="store", dest="z", type=int, default=1, help="Simulation number")
    parser.add_argument('-b', action="store", dest="amplitude", type=float, default=10, help="log(Amplitude)")
    parser.add_argument('-w', action="store", dest="width", type=float, default=1., help="width factor")
    parser.add_argument('-v', action="store", dest="bumpVelFactor", type=float, default=0, help="% of nominal")
    parser.add_argument('-p', action="store", dest="iniBumpPeakPos", type=int, default=90, help="Starting center (AU)")
    parser.add_argument('-i', action="store", dest="invertBump", type=int, default="1", help="Invert the gauss bump")
    arguments = parser.parse_args()
    main(arguments)