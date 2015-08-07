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

def set_parameter(para):
	"""to set all the paramters for plotting and calculating
	including fmin, fmax, tmin, tmax etc."""
	# if given by user in command line 
	if para: 
		para = adjust_para(para)
		if not para:
			return []
		# f_flag = set_flag('filter')
		# c_flag = set_flag('cut')
		# if f_flag and (len(para) >= 8):
		if len(para) == 10:
			# para = para[:8]
			# para.append(f_flag)
			# para.append(c_flag)

			# set filter and cut flags 

			# make sure tmin != 0
			# para = adjust_para(para)
			# print para
			return para 
		elif len(para) == 8: 
			# para = adjust_para(para)
			# if para: 
				# if not filtering; set fmin/fmax to xfmin/xfmax
			para.insert(5, para[2])
			para.insert(6, para[3])
				# para.append(f_flag)
				# para.append(c_flag)
			# if para[-4] == 0: 
			# 	para[-4] = 0.1 
				# print para
			return para 
		else: 
			return []

	# if paramters not given in command
	xtmin, xtmax = set_axis('time')
	xfmin, xfmax = set_axis('freq')
	f_flag = set_flag('filter')
	if f_flag: 
		fmin, fmax = set_bound('fas')
	else: 
		# else setting plot limits as fmin and fmax 
		fmin = xfmin
		fmax = xfmax

	tmin, tmax = set_bound('resp')
	if tmin == 0: 
		tmin = 0.1

	c_flag = set_flag('cut')

	para = [xtmin, xtmax, xfmin, xfmax, f_flag, fmin, fmax, tmin, tmax, c_flag]
	return para 
# end of set_parameter

def adjust_para(para):
	"""convert number in parameter list to floats
	set flags in parameter list to boolean
	make sure tmin != 0 """
	# convert to floats
	for i in range(0, len(para)):
		if i!=4 and i!=(len(para)-1): 
			try:
				para[i] = float(para[i])
			except ValueError:
				print "[ERROR]: invalid parameter."
				return []

	# set filter flag 
	if (para[4] in ['y', 'Y']) and len(para)==10:
		para[4] = True 
	elif (para[4] in ['n', 'N']) and len(para)==8:
		para[4] = False 
	else: 
		return []

	# set cut flag
	if para[-1] in ['y', 'Y']:
		para[-1] = True 
	elif para[-1] in ['n', 'N']:
		para[-1] = False 
	else:
		return []

	# correct tmin 
	if para[-3] == 0: 
		para[-3] = 0.1 
	return para 
# end of adjust_para


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

	return xfmin, xfmax 

def set_bound(btype):
	"""setting bounds for FAS/Response"""
	msg = ['', '', '']
	if btype == 'fas':
		msg = ['FAS', 'fmin', 'fmax']
	elif btype == 'resp':
		msg = ['Response', 'tmin', 'tmax']


	lower = raw_input('== Enter ' + msg[1] + ' for ' + msg[0] + ': ')
	upper = raw_input('== Enter ' + msg[2] + ' for ' + msg[0] + ': ')
	try: 
		lower = float(lower)
		upper = float(upper)
	except ValueError:
		print "[ERROR]: invalid input type: floats or integers ONLY."
		return set_bound(btype)

	while lower >= upper:
		print "[ERROR]: fmax must be greater than fmin."
		return set_bound(btype)

	return lower, upper 

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



def plot_signals(parameter, filenames, signal1, signal2):
	"""
	This function is to plot Signals with Fourier Amplitude Spectura. 
	"""
	plt.close('all')

	file1 = filenames[0]
	file2 = filenames[1]

	xtmin = parameter[0]
	xtmax = parameter[1]
	xfmin = parameter[2]
	xfmax = parameter[3]
	f_flag = parameter[4]
	fmin = parameter[5]
	fmax = parameter[6]
	tmin = parameter[7]
	tmax = parameter[8]

	c_flag = parameter[9]

	title = 'Acceleration ONLY: '

	samples1 = signal1.samples
	samples2 = signal2.samples
	data1 = signal1.data 
	data2 = signal2.data 
	dt1 = signal1.dt 
	dt2 = signal2.dt 

	# filtering data
	if f_flag:
		data1 = s_filter(data1, dt1, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
		data2 = s_filter(data2, dt2, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 

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

	# setting period for Response
	period = get_period(tmin, tmax)
	# rsp1 = []
	# rsp2 = []

	if not c_flag:
		# if user chooses not to cut; use origina/filted data for FAS and Response 
		freq1, fas1 = FAS(data1, dt1, points, xfmin, xfmax, 3)
		freq2, fas2 = FAS(data2, dt2, points, xfmin, xfmax, 3)
		rsp1, rsp2 = cal_acc_response(period, data1, data2, dt1, dt2)

		# for p in period:
		# 	rsp1.append(max_osc_response(data1, dt1, 0.05, p, 0, 0)[-1])
		# 	rsp2.append(max_osc_response(data2, dt2, 0.05, p, 0, 0)[-1])

	else: 
		# else uses cutted data for FAS and Response 
		freq1, fas1 = FAS(c_data1, dt1, points, xfmin, xfmax, 3)
		freq2, fas2 = FAS(c_data2, dt2, points, xfmin, xfmax, 3)
		rsp1, rsp2 = cal_acc_response(period, c_data1, c_data2, dt1, dt2)

		# for p in period:
		# 	rsp1.append(max_osc_response(c_data1, dt1, 0.05, p, 0, 0)[-1])
		# 	rsp2.append(max_osc_response(c_data2, dt2, 0.05, p, 0, 0)[-1])

	# i1 = np.where(freq1>=[xfmin])[0][0]
	# i2 = np.where(freq2>=[xfmax])[0][0]


	f, axarr = plt.subplots(nrows = 1, ncols = 3, figsize = (12, 3))

	axarr[0] = plt.subplot2grid((1, 4), (0, 0), colspan=2)
	axarr[0].set_title(title)
	axarr[0].plot(t1,c_data1,'r',t2,c_data2,'b')

	plt.legend([file1, file2])
	plt.xlim(xtmin, xtmax)

	axarr[1] = plt.subplot2grid((1, 4), (0, 2), colspan=1)
	axarr[1].set_title('Fourier Amplitude Spectra') 
	# axarr[1].plot(freq1[i1:i2],fas1[i1:i2],'r',freq2[i1:i2],fas2[i1:i2],'b')
	axarr[1].plot(freq1,fas1,'r',freq2,fas2,'b')
		
	tmp_xfmin = 0 
	if xfmin < 0.5:
		tmp_xfmin = 0
	else:
		tmp_xfmin = xfmin 
	plt.xlim(tmp_xfmin, xfmax)

	axarr[2] = plt.subplot2grid((1, 4), (0, 3))
	axarr[2].set_title("Response Spectra")
	axarr[2].set_xscale('log')
	axarr[2].plot(period,rsp1,'r',period,rsp2,'b')

	plt.xlim(tmin, tmax)

	f.tight_layout()
	plt.show()
# end of plot_signals



def plot_stations(parameter, filenames, station1, station2):
	"""
	This function is to plot two lists of psignals with Fourier Amplitude Spectra. 
	station = a list of 3 psignals for three orientation. 
	"""
	dtype = ['Displacement', 'Velocity', 'Acceleration']
	orientation = ['N/S', 'E/W', 'Up/Down']

	dt1 = station1[0].dt 
	dt2 = station2[0].dt 

	file1 = filenames[0]
	file2 = filenames[1]

	xtmin = parameter[0]
	xtmax = parameter[1]
	xfmin = parameter[2]
	xfmax = parameter[3]
	f_flag = parameter[4]
	fmin = parameter[5]
	fmax = parameter[6]
	tmin = parameter[7]
	tmax = parameter[8]
	c_flag = parameter[9]

	min_i1 = int(xtmin/dt1) 
	max_i1 = int(xtmax/dt1)

	min_i2 = int(xtmin/dt2) 
	max_i2 = int(xtmax/dt2)

	# calculating Response Spectra
	rsp1 = []
	rsp2 = []
	period = get_period(tmin, tmax)

	for i in range(0, 3):
		signal1 = station1[i]
		signal2 = station2[i]
		acc_rsp1 = []
		acc_rsp2 = []
		vel_rsp1 = []
		vel_rsp2 = []
		dis_rsp1 = []
		dis_rsp2 = []
		if c_flag:
			acc1 = signal1.accel[min_i1:max_i1]
			acc2 = signal2.accel[min_i2:max_i2]
		else: 
			acc1 = signal1.accel
			acc2 = signal2.accel

		for p in period: 
			rsp = max_osc_response(acc1, dt1, 0.05, p, 0, 0)
			dis_rsp1.append(rsp[0])
			vel_rsp1.append(rsp[1])
			acc_rsp1.append(rsp[2])

			rsp = max_osc_response(acc2, dt2, 0.05, p, 0, 0)
			dis_rsp2.append(rsp[0])
			vel_rsp2.append(rsp[1])
			acc_rsp2.append(rsp[2])
		

		rsp1.append([dis_rsp1, vel_rsp1, acc_rsp1])
		rsp2.append([dis_rsp2, vel_rsp2, acc_rsp2])

	# from displacement to velocity to acceleration
	for i in range(0, 3):
		f, axarr = plt.subplots(nrows = 3, ncols = 3, figsize = (12, 9))

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

			# filtering data
			if f_flag:
				data1 = s_filter(data1, dt1, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
				data2 = s_filter(data2, dt2, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 

			# cutting signal by bounds
			c_data1 = data1[min_i1:max_i1]
			c_data2 = data2[min_i2:max_i2]

			t1 = np.arange(xtmin, xtmax, dt1)
			t2 = np.arange(xtmin, xtmax, dt2)

			points = get_points(samples1, samples2)

			if not c_flag:
				# if user chooses not to cut; use original data to calculate FAS and Response
				freq1, fas1 = FAS(data1, dt1, points, xfmin, xfmax, 3)
				freq2, fas2 = FAS(data2, dt2, points, xfmin, xfmax, 3)

			else: 
				# use cutted signals to calculate FAS
				freq1, fas1 = FAS(c_data1, dt1, points, xfmin, xfmax, 3)
				freq2, fas2 = FAS(c_data2, dt2, points, xfmin, xfmax, 3)

			axarr[j][0] = plt.subplot2grid((3, 4), (j, 0), colspan=2, rowspan=1)
			axarr[j][0].set_title(title)
			axarr[j][0].plot(t1,c_data1,'r',t2,c_data2,'b')

			if j == 0: 
				plt.legend([file1, file2])
			plt.xlim(xtmin, xtmax)

			axarr[j][1] = plt.subplot2grid((3, 4), (j, 2), rowspan=1, colspan=1)
			axarr[j][1].set_title('Fourier Amplitude Spectra') 
			axarr[j][1].plot(freq1,fas1,'r',freq2,fas2,'b')

			tmp_xfmin = 0 
			if xfmin < 0.5:
				tmp_xfmin = 0
			else:
				tmp_xfmin = xfmin 
			plt.xlim(tmp_xfmin, xfmax)

			axarr[j][2] = plt.subplot2grid((3, 4), (j, 3), rowspan=1, colspan=1)
			axarr[j][2].set_title("Response Spectra")
			axarr[j][2].set_xscale('log')
			axarr[j][2].plot(period,rsp1[j][i],'r',period,rsp2[j][i],'b')

			plt.xlim(tmin, tmax)

		f.tight_layout()

		plt.show()
# end of plot_stations

def simple_plot(parameter, filenames, station1, station2):
	"""plotting velocity for data and FAS only
	acceleration for Response"""
	dtype = ['Displacement', 'Velocity', 'Acceleration']
	orientation = ['N/S', 'E/W', 'Up/Down']

	dt1 = station1[0].dt 
	dt2 = station2[0].dt 

	file1 = filenames[0]
	file2 = filenames[1]

	xtmin = parameter[0]
	xtmax = parameter[1]
	xfmin = parameter[2]
	xfmax = parameter[3]
	f_flag = parameter[4]
	fmin = parameter[5]
	fmax = parameter[6]
	tmin = parameter[7]
	tmax = parameter[8]
	c_flag = parameter[9]

	min_i1 = int(xtmin/dt1) 
	max_i1 = int(xtmax/dt1)

	min_i2 = int(xtmin/dt2) 
	max_i2 = int(xtmax/dt2)

	period = get_period(tmin, tmax)

	f, axarr = plt.subplots(nrows = 3, ncols = 3, figsize = (12, 9))
	for i in range(0, 3):
		# title = file1 + ' ' + file2
		title = orientation[i]

		signal1 = station1[i]
		signal2 = station2[i]

		samples1 = signal1.samples
		samples2 =  signal2.samples 

		vel1 = signal1.velo
		vel2 = signal2.velo

		acc1 = signal1.accel
		acc2 = signal2.accel

		# filtering data
		if f_flag:
			vel1 = s_filter(vel1, dt1, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
			vel2 = s_filter(vel2, dt2, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
			acc1 = s_filter(acc1, dt1, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
			acc2 = s_filter(acc2, dt2, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 

		# cutting signal by bounds
		c_vel1 = vel1[min_i1:max_i1]
		c_vel2 = vel2[min_i2:max_i2]

		c_acc1 = acc1[min_i1:max_i1]
		c_acc2 = acc2[min_i2:max_i2]

		t1 = np.arange(xtmin, xtmax, dt1)
		t2 = np.arange(xtmin, xtmax, dt2)

		points = get_points(samples1, samples2)

		if c_flag: 
			freq1, fas1 =  FAS(c_vel1, dt1, points, xfmin, xfmax, 3)
			freq2, fas2 =  FAS(c_vel2, dt2, points, xfmin, xfmax, 3)

			rsp1, rsp2 = cal_acc_response(period, c_acc1, c_acc2, dt1, dt2)
		else: 
			freq1, fas1 =  FAS(vel1, dt1, points, xfmin, xfmax, 3)
			freq2, fas2 =  FAS(vel2, dt2, points, xfmin, xfmax, 3)

			rsp1, rsp2 = cal_acc_response(period, acc1, acc2, dt1, dt2)


		axarr[i][0] = plt.subplot2grid((3, 4), (i, 0), colspan=2, rowspan=1)
		axarr[i][0].set_title(title)
		axarr[i][0].plot(t1,c_vel1,'r',t2,c_vel2,'b')

		if i == 0: 
			plt.legend([file1, file2])
			plt.xlim(xtmin, xtmax)

		axarr[i][1] = plt.subplot2grid((3, 4), (i, 2), rowspan=1, colspan=1)
		axarr[i][1].set_title('Fourier Amplitude Spectra') 
		axarr[i][1].plot(freq1,fas1,'r',freq2,fas2,'b')

		tmp_xfmin = 0
		if xfmin < 0.5:
			tmp_xfmin = 0
		else:
			tmp_xfmin = xfmin 
		plt.xlim(tmp_xfmin, xfmax)

		axarr[i][2] = plt.subplot2grid((3, 4), (i, 3), rowspan=1, colspan=1)
		axarr[i][2].set_title("Response Spectra")
		axarr[i][2].set_xscale('log')
		axarr[i][2].plot(period,rsp1,'r',period,rsp2,'b')

		plt.xlim(tmin, tmax)
	f.tight_layout()
	plt.show()
# end of simple_plot


# def test(psignal):
# 	"""to test with psignal's data"""
# 	print psignal.data[:,0], psignal.data[:,1], psignal.data[:,2]
# 	print psignal.accel
# 	print psignal.velo
# 	print psignal.displ
# 	print (psignal.data[:,0]==psignal.displ).all()
# 	print (psignal.data[:,1]==psignal.velo).all()
# 	print (psignal.data[:,2]==psignal.accel).all()
# 	print psignal.data.shape[1]


def compare_txt(parameter, file1, file2):
	if not parameter:
		print "[ERROR]: error with parameters for ploting and computing."
		return 

	signal1 = read_txt(file1)
	signal2 = read_txt(file2)

	if (not isinstance(signal1, seism_signal)) or (not isinstance(signal2, seism_signal)):
		print "[ERROR]: Invalid instance type: can only compare signal objects."
		return 

	filenames = [file1, file2]
	plot_signals(parameter, filenames, signal1, signal2)
# end of compare_txt



def compare_her(parameter, file1, file2, s_flag):
	if not parameter:
		print "[ERROR]: error with parameters for ploting and computing."
		return 

	# station = [psignal_ns, psignal_ew, psignal_up]
	station1 = read_her(file1) 
	station2 = read_her(file2)
	if len(station1) != 3 or len(station2) != 3:
		print "[ERROR]: plot_stations only handles stations with 3 channels."
		return 

	filenames = [file1, file2]

	# calling simple compare
	if s_flag:
		simple_plot(parameter, filenames, station1, station2)
	else: 
		plot_stations(parameter, filenames, station1, station2)
# end of compare_her



