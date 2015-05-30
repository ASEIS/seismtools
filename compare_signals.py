#!/usr/bin/env python
# ==========================================================================
# The program is to read two files and plot their data in a single graph. 
# It supports both .txt file and .her file. 
# 
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename): 
	"""
	The function is to read general 1-column text files. 
	"""
	f = open('SampleOutputs/' + filename, 'r')
	samples = 0 
	dt = 0.0 
	data = []
	for line in f:
		# get header 
		if "#" in line: 
			tmp = line.split()
			samples = int(tmp[4])
			dt = float(tmp[6])
		# get data 
		else:
			# print line.split()[0]
			data.append(float(line))
	data = np.array(data)
	f.close()
	return samples, dt, data 

def read_her_file(filename):
	"""
	The function is to read 10-column .her files. 
	"""
	time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt('SampleOutputs/' + filename, skiprows = 1, unpack = True)
	samples = dis_ns.size 
	dt = time[1]
	# TODO: return list of records or station object 
	return samples, dt, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up



def plot(title, samples1, dt1, data1, samples2, dt2, data2):
	plt.title(title)
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)
	plt.plot(t1,data1,'r',t2,data2,'b')
	plt.show()

# samples1, dt1, data1 = read_file("NC.NHC.V1E.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2E.txt")
# plot('Acceleration in E/W', samples1, dt1, data1, samples2, dt2, data2)

# samples1, dt1, data1 = read_file("NC.NHC.V1N.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2N.txt")
# plot('Acceleration in N/S', samples1, dt1, data1, samples2, dt2, data2)

# samples1, dt1, data1 = read_file("NC.NHC.V1Z.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2Z.txt")
# plot('Acceleration in Up/Down', samples1, dt1, data1, samples2, dt2, data2)

samples1, dt1, dis_ns1, dis_ew1, dis_up1, vel_ns1, vel_ew1, vel_up1, acc_ns1, acc_ew1, acc_up1 = read_her_file("NC.NHC.V1.her")
samples2, dt2, dis_ns2, dis_ew2, dis_up2, vel_ns2, vel_ew2, vel_up2, acc_ns2, acc_ew2, acc_up2 = read_her_file("NC.NHC.V2.her")
plot('Displacement in N/S (HER)', samples1, dt1, dis_ns1, samples2, dt2, dis_ns2) #displacement 
plot('Displacement in E/W (HER)', samples1, dt1, dis_ew1, samples2, dt2, dis_ew2) 
plot('Displacement in Up/Down (HER)', samples1, dt1, dis_up1, samples2, dt2, dis_up2) 

plot('Velocity in N/S (HER)', samples1, dt1, vel_ns1, samples2, dt2, vel_ns2) #velocity 
plot('Velocity in E/W (HER)', samples1, dt1, vel_ew1, samples2, dt2, vel_ew2) 
plot('Velocity in Up/Down (HER)', samples1, dt1, vel_up1, samples2, dt2, vel_up2) 

plot('Acceleration in N/S (HER)', samples1, dt1, acc_ns1, samples2, dt2, acc_ns2) #acceleration  
plot('Acceleration in E/W (HER)', samples1, dt1, acc_ew1, samples2, dt2, acc_ew2) 
plot('Acceleration in Up/Down (HER)', samples1, dt1, acc_up1, samples2, dt2, acc_up2) 

