from PIL import Image
from os import getcwd
import argparse

localDir = getcwd()


def main(args):
    i = 184
    rows = 6
    cols = 3
    dataDir = ""
    if args.z != 0:
        dataDir = '/simplots/'
        output = 'same' + str(args.z)
        stitchitsame(args.z, dataDir, 3, 1).save(localDir + '/simplots/collages/' + output + '.png')
        return
    if args.random:
        rows = 4
        cols = 1
        dataDir = '/simplots/sdr/r'
        output = 'random'
    if args.sdr:
        dataDir = '/simplots/sdr/r'
        output = 'sdr'
    elif args.mass:
        dataDir = '/simplots/mass/m'
        output = 'mass'
    elif args.dist:
        dataDir = '/simplots/dist/d'
        output = 'dist'
    if args.moving:
        i = 202
        cols = 4
        rows = 4
        output += '_moving'
    
    stitchit(i, dataDir, rows, cols).save(localDir + '/simplots/collages/' + output + '.png')
    

def stitchit(i, dataDir, row, column):
    im = Image.open(localDir + dataDir + str(i) +'.png')
    dst = Image.new('RGB', (im.width * column, im.height * row))
    pic = [208, 232,233, 234]
    i = 0
    for r in range(row):
        for c in range(column):
            im = Image.open(localDir + dataDir + str(pic[i]) +'.png')
            dst.paste(im, (c * im.width, r * im.height))
            i = i+1
    return dst


def stitchitsame(i, dataDir, row=3, column=1):
    ims = Image.open(localDir + dataDir + 'sdr/r' + str(i) +'.png')
    imm = Image.open(localDir + dataDir + 'mass/m' + str(i) + '.png')
    imd = Image.open(localDir + dataDir + 'dist/d' + str(i) + '.png')
    dst = Image.new('RGB', (ims.width * column, ims.height * row))

    dst.paste(imm, (0 * ims.width, 0 * ims.height))
    dst.paste(ims, (0 * ims.width, 1 * ims.height))
    dst.paste(imd, (0 * ims.width, 2 * ims.height))
    return dst


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action="store", dest="sdr", type=int, default=0, help="SDR")
    parser.add_argument('-m', action="store", dest="mass", type=int, default=0, help="Mass")
    parser.add_argument('-d', action="store", dest="dist", type=int, default=0, help="Distribution")
    parser.add_argument('-v', action="store", dest="moving", type=int, default=0, help="Moving bump")
    parser.add_argument('-x', action="store", dest="random", type=int, default=0, help="Moving bump")
    parser.add_argument('-z', action="store", dest="z", type=int, default=0, help="M")
    arguments = parser.parse_args()
    main(arguments)
