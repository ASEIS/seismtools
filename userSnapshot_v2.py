#!/usr/bin/env python
# ==========================================================================
# The program is to load data stored in plantedisplacement file, then to plot
# figures based on user's choises.
# version: August 11, 2015.
# ==========================================================================
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from userInput import Input
from planeData import *
from htools import show_progress

def components(magSelect, dataX, dataY, dataZ):
	"""decide components to use for displacement plotting"""
	magDic = {'x': dataX, 'y': dataY, 'z': dataZ}
	if len(magSelect) == 1:
		magnitude = magDic[magSelect[0]]
	elif len(magSelect) == 2:
		magnitude = np.sqrt(np.power(magDic[magSelect[0]], 2) + np.power(magDic[magSelect[1]], 2))
	elif len(magSelect) == 3:
		magnitude = np.sqrt(np.power(magDic[magSelect[0]], 2) + np.power(magDic[magSelect[1]], 2)
			+ np.power(magDic[magSelect[2]], 2))
	return magnitude
# end of disComponents

def cumulativePeak(peak, magnitude):
	"""return the peak value based on original peak and given magnitude"""
	peak = np.maximum(peak, magnitude)
	return peak
# end of cumulativeMag

def notCum(peak, data_mag):
	return data_mag

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

	dis = Data('d', disX, disY, disZ)
	# return disX, disY, disZ
	return dis

def zero_matrix(stepAlongStrike, alongStrike, stepDownDip, downDip):
	"""generate a matrix contains all zeros"""
	y = np.array(range(0, stepAlongStrike*alongStrike, stepAlongStrike))
	x = np.array(range(0, stepDownDip*downDip, stepDownDip))
	x, y = np.meshgrid(x, y)
	zeros = np.zeros_like(x)
	return zeros
# end of init_peak

def plot(peak, userInput, index):
	if userInput.barChoice == True:
		im = plt.imshow(peak, vmin=userInput.barMin,
			vmax=userInput.barMax, cmap=userInput.colorMap)
	else:
		im = plt.imshow(peak, cmap=userInput.colorMap)
	# im = plt.imshow(peak)

	plt.axis('off')
	plt.gca().invert_yaxis()

	plt.colorbar(im)
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.suptitle('t = ' + (str)(index*userInput.numSnapshots), fontsize=20)
	plt.axis('scaled')
	saveImage(userInput, index)
	plt.show()
# end of plot

def saveImage(userInput, index):
	filename = gen_fileName(userInput, index, 'png')
	plt.savefig(filename)
# end of saveImage

def gen_fileName(userInput, index, fileType):
	"""define filename based on user's inputs"""
	type_dict = {'a':'acc', 'v':'vel', 'd': 'dis'}
	mag_dict = {True: '-mag', False: ''}
	cum_dict = {True: '-cum', False: ''}

	if len(userInput.magSelect) == 1:
		magSelect = userInput.magSelect[0] + mag_dict[userInput.magnitude]
	else:
		magSelect = ''.join(userInput.magSelect)
	if userInput.snapshots == 's':
		filename = type_dict[userInput.plotType] + '-' + magSelect + cum_dict[userInput.cumulative] + '-single' + '.' + fileType
	elif userInput.snapshots == 'm':
		filename = type_dict[userInput.plotType] + '-' + magSelect + cum_dict[userInput.cumulative] + '-' + str(index*userInput.numSnapshots) + 's' + '.' + fileType

	return filename
# end of gen_fileName

def saveDat(userInput, plotData, index):
	"""print the data used to plot in a separate file"""
	filename = gen_fileName(userInput, index, 'dat')

	try:
		f = open(filename, 'w')
	except IOError, e:
		print e

	descriptor = '{:>12}'*2 + '{:>12.7f}' + '\n'
	x = np.arange(0, userInput.downDip*userInput.stepDownDip, userInput.stepDownDip, dtype=np.int)
	for i in range(0, len(plotData)):
		y = np.empty(len(plotData[i]), dtype = np.int)
		y.fill(i*userInput.stepAlongStrike)
		values = plotData[i]
		for c0, c1, c2 in zip(x, y, values):
			f.write(descriptor.format(c0, c1, c2))
	f.close()
# end of saveDat

def notSaveDat(userInput, plotData, index):
	pass

def signed(dis):
	return dis

def unsigned(dis):
	"""take the absolute value of given data"""
	dis.dataX =  np.absolute(dis.dataX)
	dis.dataY =  np.absolute(dis.dataY)
	dis.dataZ =  np.absolute(dis.dataZ)
	return dis

def derivative(data0, data, dt):
	dataX = (data.dataX - data0.dataX)/dt
	dataY = (data.dataY - data0.dataY)/dt
	dataZ = (data.dataZ - data0.dataZ)/dt
	newData = Data(data.dtype, dataX, dataY, dataZ)

	return newData
# end of derivative

def processDis(planeData, userInuput):
	dis = planeData.dis
	mag = components(userInuput.magSelect, dis.dataX, dis.dataY, dis.dataZ)
	return planeData, mag

def processVel(planeData, userInuput):
	vel = derivative(planeData.pre_dis, planeData.dis, userInuput.deltaT)
	planeData.pre_dis = planeData.dis # update planeData
	mag = components(userInuput.magSelect, vel.dataX, vel.dataY, vel.dataZ)
	return planeData, mag

def processAcc(planeData, userInuput):
	vel = derivative(planeData.pre_dis, planeData.dis, userInuput.deltaT)
	acc = derivative(planeData.pre_vel, vel, userInuput.deltaT)
	planeData.pre_dis = planeData.dis # update planeData
	planeData.pre_vel = vel
	mag = components(userInuput.magSelect, acc.dataX, acc.dataY, acc.dataZ)
	return planeData, mag


def userSnapshot(userInput):
	simulationTime = userInput.simulationTime
	deltaT = userInput.deltaT
	runtime = int(simulationTime/deltaT)
	plotType = userInput.plotType
	alongStrike = userInput.alongStrike
	downDip = userInput.downDip
	stepAlongStrike = userInput.stepAlongStrike
	stepDownDip = userInput.stepDownDip
	magSelect = userInput.magSelect
	snapshots = userInput.snapshots
	numSnapshots = userInput.numSnapshots
	magnitude = userInput.magnitude
	cumulative = userInput.cumulative

	# initializing data
	zeros = zero_matrix(stepAlongStrike, alongStrike, stepDownDip, downDip)
	peak = zeros
	data_mag = zeros
	dis = Data('d', zeros, zeros, zeros)
	dis0 = Data('d', zeros, zeros, zeros)
	vel0 = Data('v', zeros, zeros, zeros)
	planeData = PlaneData(dis, dis0, vel0)

	# define functions to process data
	process_dict = {}
	if magnitude:
		process_dict['mag'] = unsigned
	else:
		process_dict['mag'] = signed

	if plotType == 'd':
		process_dict['process'] = processDis
	elif plotType == 'v':
		process_dict['process'] = processVel
	elif plotType == 'a':
		process_dict['process'] = processAcc

	if cumulative:
		process_dict['cum'] = cumulativePeak
	else:
		process_dict['cum'] = notCum

	if userInput.printDat:
		process_dict['save'] = saveDat
	else:
		process_dict['save'] = notSaveDat


	plotType_dict = {}
	index = 0
	for i in range(0, runtime):
		planeData.dis = readFile(userInput.fp, downDip, alongStrike)

		planeData.dis = process_dict['mag'](planeData.dis)
		planeData, data_mag = process_dict['process'](planeData, userInput)
		peak = process_dict['cum'](peak, data_mag)


		if snapshots == 'm' and ((i*deltaT)%numSnapshots == 0):
			index += 1
			plot(peak, userInput, index)
			process_dict['save'](userInput, peak, index)

		# showing progress on terminal
		show_progress(i, runtime)


	# plotting only cumulative values
	if snapshots == 's':
		plot(peak, userInput, 0)
		process_dict['save'](userInput, peak, 0)
	sys.stdout.write('\n')
# end of userSnapshot


if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		userInput = Input(*argument)
	else:
		userInput = Input()

	userSnapshot(userInput)


