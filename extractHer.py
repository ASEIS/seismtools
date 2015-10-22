# /*
#  ============================================================================
# Loads the plane displacement data at given coordinates, generates corresponding
# velocity and acceleration data, in order to print the hercules file.
# Version: Sep 10, 2015
#  ============================================================================
#  */
from __future__ import division
import sys
import numpy as np
from htools import *
from scipy import interpolate
from userInput import *

def readFile(fp, downDip, alongStrike):
	"""read the binary file to get the X, Y, Z values and reshape each
	into a 2D matrix"""
	dis = np.fromfile(fp, np.float64, downDip*alongStrike*3)

	X = dis[::3] #take every third element starting at index 0
	Y = dis[1::3] #...starting at index 1
	Z = dis[2::3] #...starting at index 2

	disX = np.reshape(X, (downDip, alongStrike), order='F')
	disY = np.reshape(Y, (downDip, alongStrike), order='F')
	disZ = np.reshape(Z, (downDip, alongStrike), order='F')

	disX = disX.transpose()
	disY = disY.transpose()
	disZ = disZ.transpose()

	return disX, disY, disZ

def load_by_index(fp, alongStrike, downDip, x_coor, y_coor, num_layers):
	"""load plane data file by index, return the values of four points around the station"""
	dis = np.array([], float)
	values_x, values_y, values_z = np.array([], float), np.array([], float), np.array([], float)
	base = alongStrike*downDip*3*num_layers # number of data in past layers

	for i in range(0, len(x_coor)):
		x = x_coor[i]
		y = y_coor[i]

		# index = number of data in past layers + postion in current layer
		index = base + (downDip*x+y)*3
		offset = index*8 # offset measured in byte
		try:
			dis = np.memmap(fp, np.float64, 'r', offset, (3)) # load three numbers for x/y/z
		except ValueError:
			print "[ERROR]: unable to load file."
			sys.exit()

		values_x = np.append(values_x, dis[0])
		values_y = np.append(values_y, dis[1])
		values_z = np.append(values_z, dis[2])
	return values_x, values_y, values_z
# end of load_by_index

def bilinear_interp(x, y, data):
	"""perform bilinear interpolation"""
	values = np.array([], float)
	x0 = int(x)
	x1 = x0+1
	y0 = int(y)
	y1 = y0+1
	x_coor = np.array([x0, x1, x0, x1], dtype = int)
	y_coor = np.array([y0, y0, y1, y1], dtype = int)
	for i in range(0, x_coor.size):
		values = np.append(values, data[x_coor[i], y_coor[i]])
	f = interpolate.interp2d(x_coor, y_coor, values, kind='linear')
	new_value = f(x, y)
	return new_value
# end of bilinear_interp

def print_her(filename, dt, disData, velData, accData):
	try:
		f = open(filename, 'w')
	except IOError, e:
		print e
	dis_x = disData[0].tolist()
	vel_x = velData[0].tolist()
	acc_x = accData[0].tolist()
	dis_y = disData[1].tolist()
	vel_y = velData[1].tolist()
	acc_y = accData[1].tolist()
	dis_z = disData[2].tolist()
	vel_z = velData[2].tolist()
	acc_z = accData[2].tolist()


	# get a list of time incremented by dt
	time = [0.000]
	samples = disData[0].size
	tmp = samples

	while tmp > 1:
		time.append(time[len(time)-1] + dt)
		tmp -= 1

	# f.write('# '+str(dt)+'\n')

	descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
	f.write(descriptor.format("# time", "dis_x", "dis_y", "dis_z", "vel_x", "vel_y", "vel_z", "acc_x", "acc_y", "acc_z")) # header

	descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
	for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_x, dis_y, dis_z, vel_x, vel_y, vel_z, acc_x, acc_y, acc_z):
		f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
	f.close()
	print "*Generated hercules file at: " + filename
# end of print_her


if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		e = ExtractInput(*argument)
	else:
		e = ExtractInput()
	fp = e.fp
	simulationTime = e.simulationTime
	deltaT = e.deltaT
	alongStrike = e.alongStrike
	downDip = e.downDip
	stepAlongStrike = e.stepAlongStrike
	stepDownDip = e.stepDownDip
	x = e.x
	y = e.y

	index_x = x/stepAlongStrike
	index_y = y/stepDownDip

	# get the coordinates of points around station for bilinear interpolation
	x0 = int(index_x)
	x1 = x0+1
	y0 = int(index_y)
	y1 = y0+1
	x_coor = np.array([x0, x1, x0, x1], dtype = int)
	y_coor = np.array([y0, y0, y1, y1], dtype = int)


	runtime = int(simulationTime/deltaT)
	disX = np.array([],float)
	disY = np.array([],float)
	disZ = np.array([],float)
	for i in range(0, runtime):
		values_x, values_y, values_z = load_by_index(fp, alongStrike, downDip, x_coor, y_coor, i)
		fx = interpolate.interp2d(x_coor, y_coor, values_x, kind='linear')
		fy = interpolate.interp2d(x_coor, y_coor, values_y, kind='linear')
		fz = interpolate.interp2d(x_coor, y_coor, values_z, kind='linear')
		disX = np.append(disX, fx(index_x, index_y))
		disY = np.append(disY, fy(index_x, index_y))
		disZ = np.append(disZ, fz(index_x, index_y))


		# load the whole data; save for testing
		# dataX, dataY, dataZ = readFile(fp, downDip, alongStrike)

		# disX = np.append(disX, dataX[int(index_x), int(index_y)])
		# disY = np.append(disY, dataY[int(index_x), int(index_y)])
		# disZ = np.append(disZ, dataZ[int(index_x), int(index_y)])

		# disX = np.append(disX, bilinear_interp(index_x, index_y, dataX))
		# disY = np.append(disY, bilinear_interp(index_x, index_y, dataY))
		# disZ = np.append(disZ, bilinear_interp(index_x, index_y, dataZ))

		# showing progress on terminal
		show_progress(i, runtime)
	sys.stdout.write('\n')

	velX = derivative(disX, deltaT)
	velY = derivative(disY, deltaT)
	velZ = derivative(disZ, deltaT)

	accX = derivative(velX, deltaT)
	accY = derivative(velY, deltaT)
	accZ = derivative(velZ, deltaT)

	print_her(e.out_path, deltaT, [disX, disY, disZ], [velX, velY, velZ], [accX, accY, accZ])
# end of __main__
