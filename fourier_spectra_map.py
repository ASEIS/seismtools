# /*
#  ============================================================================
# Loads the plane displacement data to produce acceleration signals at each grid
# point, then generates a fourier spectra map of signals.
# Version: Oct 16, 2015.
#  ============================================================================
#  */
import numpy as np
import sys
# sys.path.insert(0, '/Users/kelicheng/seismtools') # insert the path to seismtools to import stools program.
from stools import smooth
from htools import plot, saveDat, derivative, dis_to_acc, show_progress, loadFile
from userInput import FourierInput
np.seterr(divide='ignore', invalid='ignore')
def not_process(data, deltaT):
	return data

def FAS(data, dt, points, freq, s_factor):
  afs = abs(np.fft.fft(data, points))*dt
  afs = smooth(afs, s_factor)
  deltaf = (1/dt)/points
  index = int(freq/deltaf)
  return afs[index]
# end of FAS



def fourier_spectra_map(userInput):
	compoDic = {'x':0, 'y':1, 'z':2}
	typeDic = {'a':dis_to_acc, 'v':derivative, 'd':not_process}
	fp = userInput.fp
	dimensionX = userInput.dimensionX
	dimensionY = userInput.dimensionY
	spaceX = userInput.spaceX
	spaceY = userInput.spaceY
	block_size = userInput.size

	simulationTime = userInput.simulationTime
	deltaT = userInput.deltaT
	frequence = userInput.frequence
	component = userInput.component
	dataType = userInput.plotType
	colorMap = userInput.colorMap


	points = int(simulationTime/deltaT)
	numGridX = dimensionX/spaceX+1
	numGridY = dimensionY/spaceY+1


	fourier = np.empty((numGridY, numGridX))
	dis = np.empty((points, block_size))
	# iterate through each grid points
	for i in range(0, numGridY):
		j = 0
		while j < numGridX:
			if j + block_size > numGridX:
				tmp_size = numGridX-j
			else:
				tmp_size = block_size

			# load the data at each time stamp
			for k in range(0, points):
				data = loadFile(fp, numGridY, numGridX, k, i, j, tmp_size)
				dis[k] = data[compoDic[component]::3]

			# for each signal, generate acceleration signal and calculate response
			for k in range(0, tmp_size):
				data = typeDic[dataType](dis[:,k], deltaT) # convert dis data to specified data type
				fourier[i][j+k] = FAS(data, deltaT, points, frequence, 3)

			j+=block_size
		show_progress(i, numGridY)
	sys.stdout.write('\n')
	if userInput.printDat:
		saveDat(userInput.out_path, dimensionX, spaceX, spaceY, fourier)
	# plot(response, colorMap)
	plot(fourier, userInput)
# end of response_spectra_map2



if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		userInput = FourierInput(*argument)
	else:
		userInput = FourierInput()

	fourier_spectra_map(userInput)





