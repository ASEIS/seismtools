#!/usr/bin/env python
# ==========================================================================
# The program is to read two files and plot their data in a single graph 
# 
# 
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename): 
	f = open(filename, 'r')
	samples = 0 
	dt = 0.0 
	data = []
	for line in f:
		# get header 
		if "#" in line: 
			tmp = line.split()
			samples = int(tmp[5])
			dt = float(tmp[7])
		# get data 
		else:
			data.append(float(line))
	data = np.array(data)
	f.close()
	return samples, dt, data 


def plot(samples1, dt1, data1, samples2, dt2, data2):
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)
	plt.plot(t1,data1,'r',t2,data2,'b')
	plt.show()

samples1, dt1, data1 = read_file("NC.NHC.V1E.txt")
samples2, dt2, data2 = read_file("NC.NHC.V2E.txt")
plot(samples1, dt1, data1, samples2, dt2, data2)

samples1, dt1, data1 = read_file("NC.NHC.V1N.txt")
samples2, dt2, data2 = read_file("NC.NHC.V2N.txt")
plot(samples1, dt1, data1, samples2, dt2, data2)

samples1, dt1, data1 = read_file("NC.NHC.V1Z.txt")
samples2, dt2, data2 = read_file("NC.NHC.V2Z.txt")
plot(samples1, dt1, data1, samples2, dt2, data2)
