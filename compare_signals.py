#!/usr/bin/env python
# ==========================================================================
# The program is to read two files and plot their data in a single graph. 
# It supports both .txt file and .her file. 
# 
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
# import math 
from seism import *
from stools import *

def set_axis(xtype):
	"""setting bounds for ploting"""
	msg = ['', '']

	if xtype == 'freq':
		msg = ['F-min', 'F-max']
	elif xtype == 'time':
		msg = ['T-min', 'T-max']

	xfmin = raw_input('== Enter ' + msg[0] + ' for ploting: ')
	xfmax = raw_input('== Enter ' + msg[1] + ' for ploting: ')

	try: 
		xfmin = float(xfmin)
		xfmax = float(xfmax)
	except ValueError:
		print "[ERROR]: invalid input type: floats or integers ONLY."
		return set_axis(xtype)

	while xfmin >= xfmax:
		print "[ERROR]: max must be greater than min."
		return set_axis(xtype)
	if xfmin < 0.1: 
		xfmin = 0

	return xfmin, xfmax 

def set_bound():
	fmin = raw_input('== Enter fmin for FAS: ')
	fmax = raw_input('== Enter fmax for FAS: ')
	try: 
		fmin = float(fmin)
		fmax = float(fmax)
	except ValueError:
		print "[ERROR]: invalid input type: floats or integers ONLY."
		return set_bound()

	while fmin >= fmax:
		print "[ERROR]: fmax must be greater than fmin."
		return set_bound()

	return fmin, fmax 

def set_flag(ftype):
	"""get user's choice on cutting signal/filter data"""
	question = ''
	if ftype == 'filter':
		question = '== Do you want to filter the data (y/n): '
	elif ftype == 'cut':
		question = '== Do you want to cut the signal (y/n): '

	flag = ''
	while flag not in ['Y', 'N']: 
		flag = raw_input(question)
		if not flag:
			return set_flag(ftype)
		flag = flag[0].upper()
	
	if flag == 'Y':
		flag = True 
	else:
		flag = False 
	return flag

# ---------------------------------------------------------------------------------------------


def read_txt(filename): 
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

def read_her(filename):
	"""
	The function is to read 10-column .her files. Return a list of psignals for each orientation/channel. 
	"""
	time = dis_ns = dis_ew = dis_up = vel_ns = vel_ew = vel_up = acc_ns = acc_ew = acc_up = np.array([],float)

	try:
		time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, comments='#', unpack = True)
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



def plot_signals(filenames, signal1, signal2):
	"""
	This function is to plot Signals with Fourier Amplitude Spectura. 
	"""
	xtmin, xtmax = set_axis('time')
	# t_flag = cut()

	xfmin, xfmax = set_axis('freq')
	fmin, fmax = set_bound()
	f_flag = set_flag('filter')
	c_flag = set_flag('cut')

	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only plot signal objects."
		return 

	plt.close('all')

	file1 = filenames[0]
	file2 = filenames[1]

	title = 'Acceleration ONLY: '

	samples1 = signal1.samples
	samples2 = signal2.samples
	data1 = signal1.data 
	data2 = signal2.data 
	dt1 = signal1.dt 
	dt2 = signal2.dt 
	t1 = np.arange(0, samples1*dt1, dt1)
	t2 = np.arange(0, samples2*dt2, dt2)

	# filtering data
	if f_flag:
		data1 = bandpass_filter(data1, dt1, fmin, fmax)
		data2 = bandpass_filter(data2, dt2, fmin, fmax)

	# cutting signals by bounds
	min_i = int(xtmin/dt1) 
	max_i = int(xtmax/dt1)
	c_data1 = data1[min_i:max_i]

	min_i = int(xtmin/dt2) 
	max_i = int(xtmax/dt2)
	c_data2 = data2[min_i:max_i]

	t1 = np.arange(xtmin, xtmax, dt1)
	t2 = np.arange(xtmin, xtmax, dt2)

	points = get_points(samples1, samples2)
	if not c_flag:
		# if user chooses not to cut; use origina/filted data for FAS 
		freq1, fas1 = FAS(data1, dt1, points, fmin, fmax, 3)
		freq2, fas2 = FAS(data2, dt2, points, fmin, fmax, 3)

	else: 
		# else uses cutted data for FAS
		freq1, fas1 = FAS(c_data1, dt1, points, fmin, fmax, 3)
		freq2, fas2 = FAS(c_data2, dt2, points, fmin, fmax, 3)

	# gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
	f, axarr = plt.subplots(nrows = 1, ncols = 2, figsize = (9, 3))

	axarr[0] = plt.subplot2grid((1, 3), (0, 0), colspan=2)
	axarr[0].set_title(title)
	axarr[0].plot(t1,c_data1,'r',t2,c_data2,'b')

	plt.legend([file1, file2])
	plt.xlim(xtmin, xtmax)

	axarr[1] = plt.subplot2grid((1, 3), (0, 2))
	axarr[1].set_title('Fourier Amplitude Spectra') 
	axarr[1].plot(freq1,fas1,'r',freq2,fas2,'b')

	plt.xlim(xfmin, xfmax)

	f.tight_layout()
	plt.show()
# end of plot_signals



def plot_stations(filenames, station1, station2):
	"""
	This function is to plot two lists of psignals with Fourier Amplitude Spectra. 
	station = a list of 3 psignals for three orientation. 
	"""
	xtmin, xtmax = set_axis('time')
	xfmin, xfmax = set_axis('freq')
	fmin, fmax = set_bound()
	f_flag = set_flag('filter')
	c_flag = set_flag('cut')

	dtype = ['Displacement', 'Velocity', 'Acceleration']
	orientation = ['N/S', 'E/W', 'Up/Down']

	if len(station1) != 3 or len(station2) != 3:
		print "[ERROR]: plot_stations only handles stations with 3 channels."
		return 

	file1 = filenames[0]
	file2 = filenames[1]

	# from displacement to velocity to acceleration
	for i in range(0, 3):
		f, axarr = plt.subplots(nrows = 3, ncols = 2, figsize = (9, 9))
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

			# filtering data
			if f_flag:
				data1 = bandpass_filter(data1, dt1, fmin, fmax)
				data2 = bandpass_filter(data2, dt2, fmin, fmax)

			# cutting signal by bounds
			min_i = int(xtmin/dt1) 
			max_i = int(xtmax/dt1)
			c_data1 = data1[min_i:max_i]

			min_i = int(xtmin/dt2) 
			max_i = int(xtmax/dt2)
			c_data2 = data2[min_i:max_i]

			t1 = np.arange(xtmin, xtmax, dt1)
			t2 = np.arange(xtmin, xtmax, dt2)

			points = get_points(samples1, samples2)
			if not c_flag:
				# if user chooses not to cut; use original data to calculate FAS
				freq1, fas1 = FAS(data1, dt1, points, fmin, fmax, 3)
				freq2, fas2 = FAS(data2, dt2, points, fmin, fmax, 3)
			else: 
				# use cutted signals to calculate FAS
				freq1, fas1 = FAS(c_data1, dt1, points, fmin, fmax, 3)
				freq2, fas2 = FAS(c_data2, dt2, points, fmin, fmax, 3)

			# 	t1 = np.arange(0, samples1*dt1, dt1)
			# 	t2 = np.arange(0, samples2*dt2, dt2)

			# axarr[j][0].set_title(title)
			# axarr[j][0].plot(t1,data1,'r',t2,data2,'b')
			# axarr[j][1].set_title('Fourier Amplitude Spectra')
			# axarr[j][1].plot(freq1,fas1,'r',freq2,fas2,'b')


			axarr[j][0] = plt.subplot2grid((3, 3), (j, 0), colspan=2, rowspan=1)
			axarr[j][0].set_title(title)
			axarr[j][0].plot(t1,c_data1,'r',t2,c_data2,'b')

			plt.legend([file1, file2])
			plt.xlim(xtmin, xtmax)

			axarr[j][1] = plt.subplot2grid((3, 3), (j, 2), rowspan=1)
			axarr[j][1].set_title('Fourier Amplitude Spectra') 
			axarr[j][1].plot(freq1,fas1,'r',freq2,fas2,'b')

			plt.xlim(xfmin, xfmax)

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
	signal1 = read_txt(file1)
	signal2 = read_txt(file2)
	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only compare signal objects."
		return 

	filenames = [file1, file2]
	plot_signals(filenames, signal1, signal2)
# end of compare_txt



def compare_her(file1, file2):
	# station = [psignal_ns, psignal_ew, psignal_up]
	station1 = read_her(file1) 
	station2 = read_her(file2)

	filenames = [file1, file2]
	plot_stations(filenames, station1, station2)

	# test(station1[0])
# end of compare_her

