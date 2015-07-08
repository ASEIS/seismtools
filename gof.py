#!/usr/bin/env python
# ===================================================================================
# The program is to read two .her files; process their signals; 
# calculate their scores with different sample rates; 
# and generate 3D matrix for scores. 
# ===================================================================================
from __future__ import division
import sys
import os
import numpy as np
import math 
from seism import *
from stools import *


destination = ''
f0 = 0.05
f1 = 0.1
f2 = 0.25 
f3 = 0.5
f4 = 1
f5 = 2
f6 = 4 
bands = [f0, f1, f2, f3, f4, f5, f6]


def get_files(): 
	file1 = ''
	file2 = ''
	filelist = ''
	global destination

	if len(sys.argv) == 2:
		filelist =  sys.argv[1]

	elif len(sys.argv) == 3: 
		file1 = sys.argv[1]
		file2 = sys.argv[2]

	while not file1:
		file1 = raw_input('== Enter the path of file1: ')

	while not file2:
		file2 = raw_input('== Enter the path of file2: ')


	# get the destination saving outputs 
	while not destination: 
		destination = raw_input('== Enter name of the directory to store outputs: ')

	# check existence of target directory 
	if not os.path.exists(destination):
		os.makedirs(destination)

	get_bands()

	return file1, file2
# end of get_files

def get_bands():
	"""
	The function is to allow user specify sample rates. 
	Without user input, sample rates are setting to default values. 
	"""
	global bands
	freq = []
	flag = True 

	while flag: 
		flag = False 
		freq = raw_input('== Enter the sequence of sample rates: ').replace(',', ' ').split()
		if not freq:
			#setting to default values 
			return

		if len(freq) == 1:
			flag = True 
		else: 
			bands = []
			for f in freq:
				try: 
					bands.append(float(f))
				except ValueError:
					print "[ERROR]: invalid sample rates"
					flag = True
					break 

			for i in range(0, len(bands)-1):
				if bands[i] >= bands[i+1]:
					print "[ERROR]: invalid sequence of sample rates"
					flag = True 
					break 

# enf of get_bands

def read_file(filename):
	"""
	The function is to read 10-column .her files. 
	Return a list of psignals for each orientation. 
	"""
	time = dis_ns = dis_ew = dis_up = vel_ns = vel_ew = vel_up = acc_ns = acc_ew = acc_up = np.array([],float)

	try:
		time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, skiprows = 1, unpack = True)
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
	return station 
# end of read_file

# ===========================================================================================
# def filter(psignal):
# 	"""The function is to call ellip_filter on each psignal """
# 	if not isinstance(psignal, seism_psignal):
# 		print "[ERROR]: encounter error filting psignal."
# 		return False 

# 	# order of unlabeled arguments is: (data, dt, order N, rp, rs, Wn) 
# 	psignal.accel = highpass_filter(psignal.accel, psignal.dt) 
# 	psignal.velo = highpass_filter(psignal.velo, psignal.dt) 
# 	psignal.displ = highpass_filter(psignal.displ, psignal.dt) 

# 	psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

# 	# psignal.print_attr()

# 	return psignal 

def filter_data(psignal, fmin, fmax):
	if not isinstance(psignal, seism_psignal):
		print "[ERROR]: encounter error filting psignal."
		return False 
	dt = psignal.dt 
	psignal.accel = bandpass_filter(psignal.accel, dt, fmin, fmax)
	psignal.velo = bandpass_filter(psignal.velo, dt, fmin, fmax)
	psignal.displ = bandpass_filter(psignal.displ, dt, fmin, fmax)

	psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

	return psignal



def S(p1, p2):
	# S(p1, p2) = 10*exp{-[(p1-p2)/min(p1, p2)]^2}
	# if p1 == p2:
	# 	return 10 
	s = 10*np.exp(-((p1-p2)/min(p1, p2))**2)
	return s 


def cal_peak(data1, data2):
	"""
	calculate the socres for peak acc/vel/dis.
	score = S(max|data1|, max|data2|) 
	"""
	p1 = np.amax(np.absolute(data1))
	p2 = np.amax(np.absolute(data2))
	score = S(p1, p2)
	return score 

def I(data, dt):
	# I(t) = max|integral(data^2)dt|
	return np.amax(np.cumsum(data*data)*dt)


def cal_SI(signal1, signal2):
	"""
	score = S(IA1, IA2) for Arias intensity 
	score = S(IE1, IE2) for Energy integral 
	"""
	I1 = I(signal1.accel, signal1.dt)
	I2 = I(signal2.accel, signal2.dt)
	SIa = S(I1, I2)

	I1 = I(signal1.velo, signal1.dt)
	I2 = I(signal2.velo, signal2.dt)
	SIv = S(I1, I2)
	return SIa, SIv 

def N(data, dt):
	"""N = Ie(t)/IE = Ia(t)/IA"""
	return np.cumsum(data*data)*dt/I(data, dt)

def F(N1, N2):
	return np.absolute(N1-N2)


def cal_SD(signal1, signal2):
	"""SD = 10*(1-max(F))"""
	N1 = N(signal1.accel, signal1.dt)
	N2 = N(signal2.accel, signal2.dt)
	SDa = 10*(1-np.amax(F(N1, N2)))

	N1 = N(signal1.velo, signal1.dt)
	N2 = N(signal2.velo, signal2.dt)
	SDe = 10*(1-np.amax(F(N1, N2)))

	return SDa, SDe

def cal_Sfs(signal1, signal2, fmin, fmax):
	"""
	calculate the score for Fourier Spectra
	Sfs = mean(S(FS1, FS2))
	"""
	points = get_points(signal1.samples, signal2.samples)
	fs1 = FAS(signal1.velo, signal1.dt, points, fmin, fmax, 3)[-1]
	fs2 = FAS(signal2.velo, signal2.dt, points, fmin, fmax, 3)[-1]
	s = np.array([],float)

	for i in range(0, fs1.size):
		s = np.append(s, S(fs1[i], fs2[i]))

	# print s.size 
	Sfs = np.mean(s)
	return Sfs

def cal_C(a1, a2, dt):
	"""
	calculate the score for Cross Correlation 
	C* = 10*max(C(a1, a2), 0)
	C = integral(a1, a2)dt/((integral(a1^2)dt^1/2)*(integral(a2^2)dt^1/2))
	"""
	x = np.cumsum(a1*a2)*dt
	y = np.cumsum(a1*a1)*dt
	z = np.cumsum(a2*a2)*dt
	c = x/(np.power(y, 0.5)*np.power(z, 0.5))
	cc = 10*np.amax(c, 0)
	return cc 


def cal_Ssa(signal1, signal2, fmin, fmax):
	"""Calculate the score for Reponse Spectra"""
	period = get_period(fmin, fmax)
	SA1 = []
	SA2 = []
	for p in period:
		SA1.append(max_osc_response(signal1.accel, signal1.dt, 0.05, p, 0, 0)[-1])
		SA2.append(max_osc_response(signal2.accel, signal2.dt, 0.05, p, 0, 0)[-1])

	ss = []
	for i in range(0, len(SA1)):
		ss.append(S(SA1[i], SA2[i]))

	return np.mean(ss)

def duration(signal):
	""" get the total duration of signal """
	data = signal.velo
	dt = signal.dt

	# E = max|integral(v^2)dt|
	E = I(data, dt)
	E5 = 0.05*E
	E95 = 0.95*E
	T5 = 0 
	T95 = 0 

	energy = np.cumsum(data*data)*dt
	for i in range(1, energy.size):
		if energy[i-1] <= E5 <= energy[i]:
			T5 = i 
		if energy[i-1] <= E95 <= energy[i]:
			T95 = i 
			break 

	D = (T95-T5)*dt
	return D

def cal_D(signal1, signal2):
	""" calculate the score for duration """
	D1 = duration(signal1)
	D2 = duration(signal2)
	return S(D1, D2)




# ============================================================================================

def scores_matrix(station1, station2):
	"""
	Generate the 3D matrix of scores 
	"""
	global bands 
	bands.insert(0, bands[len(bands)-1])
	# print bands 
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = c11 = avg = 0.0

	matrix = np.empty((4, len(bands)+1, 13))

	for i in range(1, len(station1)+1):
		for j in range(0, len(bands)-1): 
			if j == 0:
				# BB-Bn
				fmin = bands[j+1]
				fmax = bands[j]
			else: 
				# Bn-Bn+1
				fmin = bands[j]
				fmax = bands[j+1]
			# print fmin, fmax 


			signal1 = station1[i-1]
			signal2 = station2[i-1]

			# filtering data
			signal1 = filter_data(signal1, fmin, fmax)
			signal2 = filter_data(signal2, fmin, fmax)

			# signal1.print_attr()
			# signal2.print_attr()


			c1, c2 = cal_SD(signal1, signal2)
			c3, c4 = cal_SI(signal1, signal2)
			c5 = cal_peak(signal1.accel, signal2.accel)
			c6 = cal_peak(signal1.velo, signal2.velo)
			c7 = cal_peak(signal1.displ, signal2.displ)

			c8 = cal_Ssa(signal1, signal2, fmin, fmax)
			c9 = cal_Sfs(signal1, signal2, fmin, fmax)
			c10 = cal_C(signal1.accel, signal2.accel, signal1.dt)

			c11 = cal_D(signal1, signal2)

			scores = np.array([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11], float)
			# scores = np.append(scores, np.average(scores))
			T = ((c1+c2)/2 + (c3+c4)/2 + c5 + c6 + c7 + c8 + c9 + c10 + c11)/9 
			A = (c1+c2+c3+c4+c5+c6+c7+c8+c9+c10)/10 
			scores = np.insert(scores, 0, T)
			scores = np.insert(scores, 1, A)
			scores = np.around(scores, decimals=2)

			matrix[i][j] = scores 

		SA = np.array([],float)
		CA = np.array([],float)

		# calculate the average score of all bands 
		# FS2 = avg(b1-b6)
		# FS1 = avg(bb-b6)
		for j in range(0, 13):
			# SA = avg(B1...Bn)
			avg2 = np.average(matrix[i][:,j][1:len(bands)])
			SA = np.append(SA, avg2)
			SA = np.around(SA, decimals=2)

			# CA = avg(BB...Bn)
			avg1 = np.average(matrix[i][:,j][:len(bands)])
			CA = np.append(CA, avg1)
			CA = np.around(CA, decimals=2)


		# print fs1 
		# print fs2
		matrix[i][-1] = CA
		matrix[i][-2] = SA


	# insert the slide contain all AVERAGE values in front 
	for i in range(0, len(bands)+1):
		for j in range(0, 12):
			average = (matrix[1][i][j] + matrix[2][i][j] + matrix[3][i][j])/3
			matrix[0][i][j] = round(average, 2)
	
	return matrix

def summary(matrix):
	""" generate a summary matrix contain average scores 
		calculated with different methods. 
		including SA_A; CA_A; SA_T; CA_T """
	s = np.empty((4, 4))
	SA_A = CA_A = SA_T = CA_T = 0.0
	for i in range(0, 4):
		SA_A = matrix[i][-2][1]
		CA_A = matrix[i][-1][1]
		SA_T = matrix[i][-2][0]
		CA_T = matrix[i][-1][0]

		avg = np.array([SA_A, CA_A, SA_T, CA_T], float)
		s[i] = avg 
	return s 
# end of summary
# --------------------------------------------------------------------------------------------------


file1, file2 = get_files()
station1 = read_file(file1)
station2 = read_file(file2)
# sig1 = filt(station1[0])
# filt(station1[1])
# filt(station1[2])
# sig2 = filt(station2[0])
# filt(station2[1])
# filt(station2[2])
matrix = scores_matrix(station1, station2)


def print_scores(matrix):
	filename = "scores.txt"
	s = ''
	try:
		f = open(destination + '/' + filename, 'a')
	except IOError, e:
		print e

	for i in range(0, 3):
		for j in range(0, len(bands)):
			for k in range(0, len(matrix[i][j])):
				s += str(matrix[i][j][k]) + ' '

	f.write(s + '\n')
	f.close()
# end of print_scores

def print_matrix(path, matrix):
	""" generate the file containing the score matrix of two files. """
	# header = "# GOF " + file1 + ' ' + file2 
	s = summary(matrix)
	label = ['AVG', 'N', 'E', 'UP']
	# reading data of summary matrix by column 
	SA_A = s[:,0]
	CA_A = s[:,1]
	SA_T = s[:,2]
	CA_T = s[:,3]
	try:
		f = open(path, 'w')
	except IOError, e:
		print e
		# return 
	
	# printing summary matrix 
	descriptor = '{:>12}' + '  {:>12}'*4 + '\n'
	f.write(descriptor.format("# Total Average", "SA_A", "CA_A", "SA_T", "CA_T")) 

	descriptor = '{:>12}' + '  {:>12.2f}'*4 + '\n'
	for l, s1, c1, s2, c2 in zip(label, SA_A, CA_A, SA_T, CA_T):
		f.write(descriptor.format(l, s1, c1, s2, c2))

	f.write('# ----------------------------------------------------------------------------------------------\n')

	# generte row and column labels 
	num_b = len(matrix[0])-2 
	c_label = "BB" 
	for i in range(1, num_b):
		c_label += ',B'+str(i)
	c_label += ",SA,CA"
	c_label = c_label.split(',')
	c_label.insert(0, '')

	label1 = ['Average', 'North', 'East', 'Up']
	r_label = ['T', 'A', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']

	d1 = '{:>12}' + '  {:>12}'*(num_b+2) + '\n'
	d2 = '{:>12}' + '  {:>12.2f}'*(num_b+2)+'\n'

	# printing matrix 
	for i in range(0, len(matrix)):
		c_label[0] = '# '+label1[i]
		f.write(d1.format(*c_label)) 

		for s in zip(r_label, *matrix[i]):
			f.write(d2.format(*s))
		f.write('# ----------------------------------------------------------------------------------------------\n')
	f.close()

	pass 
# end of print_matrix


# print_scores(matrix)


for i in range(0, 4):
	for j in range(0, len(bands)+1):
		print matrix[i][j]
	print "---------------------------------------------------------------------------------------"

print_matrix("test.txt", matrix)



