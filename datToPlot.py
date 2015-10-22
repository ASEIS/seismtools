#!/usr/bin/env python
# ==================================================================
# The program reads data file(s) and generates corresponding plots.
# ==================================================================
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from userInput import *

def getMultiple(path, numSnapshots):
	"""searching for a sequence of snapshots"""
	filename = path.split('/')[-1]
	dirpath = path.replace(filename, '')
	time = filename.split('-')[-1]
	filelist = [path]

	while True:
		try:
			index = int(time.replace('s.dat', ''))
		except ValueError:
			print "[ERROR]: not a data file from multiple snapshots."
			sys.exit()
		new_time = str(index+numSnapshots) + 's.dat'

		filename = filename.replace(time, new_time)
		target = dirpath + filename

		if os.path.isfile(target):
			filelist.append(target)
			time = new_time
		else:
			break
	return filelist
# end of getMultiple

def datToMatrix(filename, alongStrike, downDip):
	try:
		x, y, data = np.loadtxt(filename, unpack = True)
	except IOError:
		print "[ERROR]: unable to laod file " + filename + '.'

	data = np.reshape(data, (downDip, alongStrike), order='F')
	data = data.transpose()
	return data


def plot(data, cmap):
	im = plt.imshow(data, cmap = cmap)
	plt.axis('off')
	plt.gca().invert_yaxis()

	plt.colorbar(im)
	plt.xlabel('X')
	plt.ylabel('Y')
	# plt.suptitle('t = ' + (str)(index), fontsize=20)
	plt.axis('scaled')
	plt.show()
# end of plot


if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		i = DatInput(*argument)
	else:
		i = DatInput()

	if i.snapshots == 's':
		data = datToMatrix(i.path, i.alongStrike, i.downDip)
		plot(data, i.colorMap)
	elif i.snapshots == 'm':
		filelist = getMultiple(i.path, i.numSnapshots)
		for f in filelist:
			data = datToMatrix(f, i.alongStrike, i.downDip)
			plot(data, i.colorMap)

