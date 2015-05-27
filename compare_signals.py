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
	f = open(filename, 'r')
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
	samples = sum(1 for line in open(filename)) - 1
	f = open(filename, 'r')
	dt = 0.0 
	dis_ns = []
	dis_ew = []
	dis_up = []
	vel_ns = []
	vel_ew = []
	vel_up = []
	acc_ns = []
	acc_ew = []
	acc_up = []
	line_num = -1 
	for line in f:
		line_num += 1 
		if line_num == 2: 
			dt = float(line.split()[0])
		if not "#" in line: 
			tmp = line.split()
			dis_ns.append(float(tmp[1]))
			dis_ew.append(float(tmp[2]))
			dis_up.append(float(tmp[3]))
			vel_ns.append(float(tmp[4]))
			vel_ew.append(float(tmp[5]))
			vel_up.append(float(tmp[6]))
			acc_ns.append(float(tmp[7]))
			acc_ew.append(float(tmp[8]))
			acc_up.append(float(tmp[9]))
	dis_ns = np.array(dis_ns)
	dis_ew = np.array(dis_ew)
	dis_up = np.array(dis_up)
	vel_ns = np.array(vel_ns)
	vel_ew = np.array(vel_ew)
	vel_up = np.array(vel_up)
	acc_ns = np.array(acc_ns)
	acc_ew = np.array(acc_ew)
	acc_up = np.array(acc_up)
	f.close()
	samples = int(line_num)
	return samples, dt, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up



def plot(samples1, dt1, data1, samples2, dt2, data2):
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)
	plt.plot(t1,data1,'r',t2,data2,'b')
	plt.show()

# samples1, dt1, data1 = read_file("NC.NHC.V1E.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2E.txt")
# plot(samples1, dt1, data1, samples2, dt2, data2)

# samples1, dt1, data1 = read_file("NC.NHC.V1N.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2N.txt")
# plot(samples1, dt1, data1, samples2, dt2, data2)

# samples1, dt1, data1 = read_file("NC.NHC.V1Z.txt")
# samples2, dt2, data2 = read_file("NC.NHC.V2Z.txt")
# plot(samples1, dt1, data1, samples2, dt2, data2)

samples1, dt1, dis_ns1, dis_ew1, dis_up1, vel_ns1, vel_ew1, vel_up1, acc_ns1, acc_ew1, acc_up1 = read_her_file("NC.NHC.V1.her")
samples2, dt2, dis_ns2, dis_ew2, dis_up2, vel_ns2, vel_ew2, vel_up2, acc_ns2, acc_ew2, acc_up2 = read_her_file("NC.NHC.V2.her")
plot(samples1, dt1, dis_ns1, samples2, dt2, dis_ns2) #displacement 
plot(samples1, dt1, dis_ew1, samples2, dt2, dis_ew2) 
plot(samples1, dt1, dis_up1, samples2, dt2, dis_up2) 

plot(samples1, dt1, vel_ns1, samples2, dt2, vel_ns2) #velocity 
plot(samples1, dt1, vel_ew1, samples2, dt2, vel_ew2) 
plot(samples1, dt1, vel_up1, samples2, dt2, vel_up2) 

plot(samples1, dt1, acc_ns1, samples2, dt2, acc_ns2) #acceleration  
plot(samples1, dt1, acc_ew1, samples2, dt2, acc_ew2) 
plot(samples1, dt1, acc_up1, samples2, dt2, acc_up2) 