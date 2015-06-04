#!/usr/bin/env python
# ==========================================================================
# The program is to read two files and plot their data in a single graph. 
# It supports both .txt file and .her file. 
# 
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt
from seism import *

def read_file(filename): 
	"""
	The function is to read general 1-column text files. Return a signal object. 
	"""
	try:
		f = open(filename, 'r')
	except IOError, e:
		print e
		return 
	
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
	return seism_signal(samples, dt, data, 'a')
	# return samples, dt, data 

def read_her_file(filename):
	"""
	The function is to read 10-column .her files. Return a list of psignals for each orientation/channel. 
	"""
	# TODO: add exception 
	time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, skiprows = 1, unpack = True)
	samples = dis_ns.size 
	dt = time[1]

	# samples, dt, data, acceleration, velocity, displacement 
	psignal_ns = seism_psignal(samples, dt, np.c_[acc_ns, vel_ns, dis_ns], 'a', acc_ns, vel_ns, dis_ns)
	psignal_ew = seism_psignal(samples, dt, np.c_[acc_ew, vel_ew, dis_ew], 'a', acc_ew, vel_ew, dis_ew)
	psignal_up = seism_psignal(samples, dt, np.c_[acc_up, vel_up, dis_up], 'a', acc_up, vel_up, dis_up)

	station = [psignal_ns, psignal_ew, psignal_up]
	# return samples, dt, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up
	return station 


# def plot(title, signal1, signal2):
# 	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
# 		print "[ERROR]: Invalid instance type: can only plot signal objects."
# 		return 
# 	plt.title(title)
# 	samples1 = signal1.samples
# 	data1 = signal1.data
# 	dt1 = signal1.dt


# 	samples2 = signal2.samples
# 	data2 = signal2.data
# 	dt2 = signal2.dt

# 	t1 = np.arange(0, samples1*dt1, dt1)
# 	t2 = np.arange(0, samples2*dt2, dt2)
# 	plt.plot(t1,data1,'r',t2,data2,'b')
# 	plt.show()



def plot(title, samples1, dt1, data1, samples2, dt2, data2):
	plt.title(title)
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)
	plt.plot(t1,data1,'r',t2,data2,'b')
	plt.show()


def compare_txt(file1, file2):
	# revert the order of files to V1, V2
	if 'V1' not in file1.split('.')[-2]: 
		tmp = file1
		file1 = file2 
		file2 = tmp 

	signal1 = read_file(file1)
	signal2 = read_file(file2)
	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only compare signal objects."
		return 
	samples1 = signal1.samples
	data1 = signal1.data
	dt1 = signal1.dt
	samples2 = signal2.samples
	data2 = signal2.data
	dt2 = signal2.dt
	plot('Acceleration: \n' + file1 + ' ' + file2, samples1, dt1, data1, samples2, dt2, data2)
# end of compare_txt


def compare_her(file1, file2):
	# revert the order of files to V1, V2
	if 'V1' not in file1.split('.')[-2]: 
		tmp = file1
		file1 = file2 
		file2 = tmp 
	# station = [psignal_ns, psignal_ew, psignal_up]
	station1 = read_her_file(file1) 
	station2 = read_her_file(file2)

	samples1 = station1[0].samples
	dt1 = station1[0].dt
	samples2 = station2[0].samples
	dt2 = station2[0].dt

	for i in range(0, 3):
		if (not isinstance(station1[i], seism_signal)) or (not isinstance(station2[i], seism_signal)):
			print "[ERROR]: Invalid instance type: can only compare psignal objects."
			return 

	plot('Displacement in N/S: \n' + file1 + ' ' + file2, samples1, dt1, station1[0].displ, samples2, dt2, station2[0].displ) #displacement 
	plot('Displacement in E/W: \n' + file1 + ' ' + file2, samples1, dt1, station1[1].displ, samples2, dt2, station2[1].displ) 
	plot('Displacement in Up/Down: \n' + file1 + ' ' + file2, samples1, dt1, station1[2].displ, samples2, dt2, station2[2].displ) 

	plot('Velocity in N/S: \n' + file1 + ' ' + file2, samples1, dt1, station1[0].velo, samples2, dt2, station2[0].velo) #velocity 
	plot('Velocity in E/W: \n' + file1 + ' ' + file2, samples1, dt1, station1[1].velo, samples2, dt2, station2[1].velo) 
	plot('Velocity in Up/Down: \n' + file1 + ' ' + file2, samples1, dt1, station1[2].velo, samples2, dt2, station2[2].velo) 

	plot('Acceleration in N/S: \n' + file1 + ' ' + file2, samples1, dt1, station1[0].accel, samples2, dt2, station2[0].accel) #acceleration  
	plot('Acceleration in E/W: \n' + file1 + ' ' + file2, samples1, dt1, station1[1].accel, samples2, dt2, station2[1].accel) 
	plot('Acceleration in Up/Down: \n' + file1 + ' ' + file2, samples1, dt1, station1[2].accel, samples2, dt2, station2[2].accel) 
# end of compare_her


	# samples1, dt1, dis_ns1, dis_ew1, dis_up1, vel_ns1, vel_ew1, vel_up1, acc_ns1, acc_ew1, acc_up1 = read_her_file(file1)
	# samples2, dt2, dis_ns2, dis_ew2, dis_up2, vel_ns2, vel_ew2, vel_up2, acc_ns2, acc_ew2, acc_up2 = read_her_file(file2)
	# plot('Displacement in N/S (HER)', samples1, dt1, dis_ns1, samples2, dt2, dis_ns2) #displacement 
	# plot('Displacement in E/W (HER)', samples1, dt1, dis_ew1, samples2, dt2, dis_ew2) 
	# plot('Displacement in Up/Down (HER)', samples1, dt1, dis_up1, samples2, dt2, dis_up2) 

	# plot('Velocity in N/S (HER)', samples1, dt1, vel_ns1, samples2, dt2, vel_ns2) #velocity 
	# plot('Velocity in E/W (HER)', samples1, dt1, vel_ew1, samples2, dt2, vel_ew2) 
	# plot('Velocity in Up/Down (HER)', samples1, dt1, vel_up1, samples2, dt2, vel_up2) 

	# plot('Acceleration in N/S (HER)', samples1, dt1, acc_ns1, samples2, dt2, acc_ns2) #acceleration  
	# plot('Acceleration in E/W (HER)', samples1, dt1, acc_ew1, samples2, dt2, acc_ew2) 
	# plot('Acceleration in Up/Down (HER)', samples1, dt1, acc_up1, samples2, dt2, acc_up2) 


# compare_her('NC.NHC.V2.her', 'NC.NHC.V1.her')


# compare_txt('NC.NHC.V1N.txt', 'NC.NHC.V2N.txt')
# compare_txt('NC.NHC.V1E.txt', 'NC.NHC.V2E.txt')
# compare_txt('NC.NHC.V1Z.txt', 'NC.NHC.V2Z.txt')