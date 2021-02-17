from PIL import Image
from os import getcwd
import argparse

localDir = getcwd()


def main(args):
	i = 184
	rows = 6
	cols = 3
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
	for r in range(row):
		for c in range(column):
			im = Image.open(localDir + dataDir + str(i) +'.png')
			dst.paste(im, (c * im.width, r * im.height))
			i = i + 1
	return dst


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action="store", dest="sdr", type=int, default=0, help="SDR")
    parser.add_argument('-m', action="store", dest="mass", type=int, default=0, help="Mass")
    parser.add_argument('-d', action="store", dest="dist", type=int, default=0, help="Distribution")
    parser.add_argument('-v', action="store", dest="moving", type=int, default=0, help="Moving bump")
    arguments = parser.parse_args()
    main(arguments)
