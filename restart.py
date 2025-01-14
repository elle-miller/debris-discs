from dustpy import readdump
from os import path
from dustpy import constants as c
import argparse
import numpy as np

# This is a script that will restart your simulation if it stops running!! eg out of space or lost connection. Eg to restart from directory 22
# $ python restart.py -z 22

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
    arguments = parser.parse_args()
    main(arguments)
