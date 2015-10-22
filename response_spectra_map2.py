# /*
#  ============================================================================
# Loads the plane displacement data to produce acceleration signals at each grid
# point, then generates a response spectra map of signals.
# Version: Oct 05, 2015.
#  ============================================================================
#  */
import numpy as np
import sys
# sys.path.insert(0, '/Users/kelicheng/seismtools') # insert the path to seismtools to import stools program.
from stools import max_osc_response
from htools import plot, saveDat, dis_to_acc, show_progress, loadFile
from userInput import ResponseInput2
np.seterr(divide='ignore', invalid='ignore')

def response_spectra_map2(userInput):
	compoDic = {'x':0, 'y':1, 'z':2}
	typeDic = {'a':2, 'v':'1', 'd':0}
	fp = userInput.fp
	dimensionX = userInput.dimensionX
	dimensionY = userInput.dimensionY
	spaceX = userInput.spaceX
	spaceY = userInput.spaceY
	block_size = userInput.size

	simulationTime = userInput.simulationTime
	deltaT = userInput.deltaT
	period = userInput.period
	component = userInput.component
	responseType = userInput.plotType
	colorMap = userInput.colorMap


	numLayer = int(simulationTime/deltaT)
	numGridX = dimensionX/spaceX+1
	numGridY = dimensionY/spaceY+1

	response = np.empty((numGridY, numGridX))
	dis = np.empty((numLayer, block_size))
	# iterate through each grid points
	for i in range(0, numGridY):
		j = 0
		while j < numGridX:
			if j + block_size > numGridX:
				tmp_size = numGridX-j
			else:
				tmp_size = block_size

			# load the data at each time stamp
			for k in range(0, numLayer):
				data = loadFile(fp, numGridY, numGridX, k, i, j, tmp_size)
				dis[k] = data[compoDic[component]::3]

			# for each signal, generate acceleration signal and calculate response
			for k in range(0, tmp_size):
				acc = dis_to_acc(dis[:,k], deltaT)
				response[i][j+k] = max_osc_response(acc, deltaT, 0.05, period, 0, 0)[typeDic[responseType]]
			j+=block_size
		show_progress(i, numGridY)
	sys.stdout.write('\n')
	if userInput.printDat:
		saveDat(userInput.out_path, dimensionX, spaceX, spaceY, response)
	# plot(response, colorMap)
	plot(response, userInput)
# end of response_spectra_map2



if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		userInput = ResponseInput2(*argument)
	else:
		userInput = ResponseInput2()

	response_spectra_map2(userInput)





