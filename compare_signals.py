#!/usr/bin/env python
# ==========================================================================
# The program is to read two files and plot their data in a single graph. 
# It supports both .txt file and .her file. 
# 
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt
# import math 
from seism import *
from stools import *

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
			samples = int(tmp[6])
			dt = float(tmp[7])
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
	time = dis_ns = dis_ew = dis_up = vel_ns = vel_ew = vel_up = acc_ns = acc_ew = acc_up = np.array([],float)

	try:
		time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, comments='#', unpack = True)
		print time
	except IOError:
		print "[ERROR]: error loading her file. "
		return False  

	samples = dis_ns.size 
	dt = time[1]

	# samples, dt, data, acceleration, velocity, displacement 
	psignal_ns = seism_psignal(samples, dt, np.c_[dis_ns, vel_ns, acc_ns], 'c', acc_ns, vel_ns, dis_ns)
	psignal_ew = seism_psignal(samples, dt, np.c_[dis_ew, vel_ew, acc_ew], 'c', acc_ew, vel_ew, dis_ew)
	psignal_up = seism_psignal(samples, dt, np.c_[dis_up, vel_up, acc_up], 'c', acc_up, vel_up, dis_up)

	station = [psignal_ns, psignal_ew, psignal_up]
	# return samples, dt, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up
	return station 



def plot_signals(title, signal1, signal2):
	"""
	This function is to plot Signals with Fourier Amplitude Spectura. 
	"""
	# global fmin
	# global fmax 
	fmin = 0.05
	fmax = 4

	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only plot signal objects."
		return 

	plt.close('all')

	samples1 = signal1.samples
	samples2 = signal2.samples
	data1 = signal1.data 
	data2 = signal2.data 
	dt1 = signal1.dt 
	dt2 = signal2.dt 
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)

	points = get_points(samples1, samples2)
	freq1, fas1 = FAS(signal1.data, dt1, points, fmin, fmax, 3)
	freq2, fas2 = FAS(signal2.data, dt2, points, fmin, fmax, 3)

	f, axarr = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5))
	axarr[0].set_title(title)
	axarr[0].plot(t1,data1,'r',t2,data2,'b')

	axarr[1].set_title('Fourier Amplitude Spectra')
	axarr[1].plot(freq1,fas1,'r',freq2,fas2,'b')
	plt.xlim(0, fmax)

	f.tight_layout()
	plt.show()
# end of plot_signals



def plot_stations(station1, station2):
	"""
	This function is to plot two lists of psignals with Fourier Amplitude Spectra. 
	station = a list of 3 psignals for three orientation. 
	"""
	# global fmin 
	# global fmax 
	fmin = 0.05
	fmax = 4

	dtype = ['Displacement', 'Velocity', 'Acceleration']
	orientation = ['N/S', 'E/W', 'Up/Down']

	if len(station1) != 3 or len(station2) != 3:
		print "[ERROR]: plot_stations only handles stations with 3 channels."
		return 

	# from displacement to velocity to acceleration
	for i in range(0, 3):
		f, axarr = plt.subplots(nrows = 3, ncols = 2, figsize = (12, 12))
		# iterative through psignals in each station 
		for j in range(0, 3):
			title = dtype[i] + ' in ' + orientation[j]
			signal1 = station1[j]
			signal2 = station2[j]
			if (not isinstance(signal1, seism_psignal)) or (not isinstance(signal2, seism_psignal)):
				print "[PLOT ERROR]: Invalid instance type: can only plot psignal objects."
				return 
			if signal1.data.shape[1] != 3 or signal2.data.shape[1] != 3:
				print "[PLOT ERROR]: psignal's data must be 3-column numpy array for displ, velo, accel."
				return 

			samples1 = signal1.samples
			samples2 = signal2.samples

			# psignal.data = [displ, velo, accel]
			data1 = signal1.data[:,i]
			data2 = signal2.data[:,i]

			dt1 = signal1.dt 
			dt2 = signal2.dt 

			t1 = np.arange(0, samples1*dt1, dt1)
			t2 = np.arange(0, samples2*dt2, dt2)

			points = get_points(samples1, samples2)
			freq1, fas1 = FAS(data1, dt1, points, fmin, fmax, 3)
			freq2, fas2 = FAS(data2, dt2, points, fmin, fmax, 3)

			axarr[j][0].set_title(title)
			axarr[j][0].plot(t1,data1,'r',t2,data2,'b')
			axarr[j][1].set_title('Fourier Amplitude Spectra')
			axarr[j][1].plot(freq1,fas1,'r',freq2,fas2,'b')
			plt.xlim(0, fmax)
			f.tight_layout()

		plt.show()
# end of plot_stations


def test(psignal):
	"""to test with psignal's data"""
	print psignal.data[:,0], psignal.data[:,1], psignal.data[:,2]
	print psignal.accel
	print psignal.velo
	print psignal.displ
	print (psignal.data[:,0]==psignal.displ).all()
	print (psignal.data[:,1]==psignal.velo).all()
	print (psignal.data[:,2]==psignal.accel).all()
	print psignal.data.shape[1]


def compare_txt(file1, file2):
	# revert the order of files to V1, V2
	if 'V1' not in file1.split('.')[-2]: 
		tmp = file1
		file1 = file2 
		file2 = tmp 

	print file1, file2

	signal1 = read_file(file1)
	signal2 = read_file(file2)
	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only compare signal objects."
		return 

	title = 'Acceleration ONLY: \n' + file1 + ' ' + file2
	plot_signals(title, signal1, signal2)
# end of compare_txt



def compare_her(file1, file2):
	# revert the order of files to V1, V2
	if 'V1' not in file1.split('.')[-2]: 
		tmp = file1
		file1 = file2 
		file2 = tmp 

	print file1, file2

	# station = [psignal_ns, psignal_ew, psignal_up]
	station1 = read_her_file(file1) 
	station2 = read_her_file(file2)

	plot_stations(station1, station2)

	# test(station1[0])
# end of compare_her

