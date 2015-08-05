#!/usr/bin/env python
# ===================================================================================
# The program receives two stations, each contains three signals,  
# then calculates signals' scores with different sample rates; 
# and generate 3D matrix for scores. 
# ===================================================================================
from __future__ import division
import sys
import os
import numpy as np
import math 
from seism import *
from stools import *

np.seterr(divide='ignore', invalid='ignore')


def update():
	"""showing progress"""
	sys.stdout.write('-')
	sys.stdout.flush()


def filter_data(psignal, fmin, fmax):
	if not isinstance(psignal, seism_psignal):
		print "[ERROR]: encounter error filting psignal."
		return False 
	dt = psignal.dt 

	psignal.accel = s_filter(psignal.accel, dt, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
	psignal.velo = s_filter(psignal.velo, dt, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
	psignal.displ = s_filter(psignal.displ, dt, type = 'bandpass', family = 'ellip', fmin = fmin, fmax = fmax, N = 3, rp = 0.1, rs = 100) 
	
	psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

	return psignal

def S(p1, p2):
	# S(p1, p2) = 10*exp{-[(p1-p2)/min(p1, p2)]^2}
	if min(p1, p2) == 0:
		return 10 
	s = 10*np.exp(-((p1-p2)/min(p1, p2))**2)
	return s 


def cal_peak(data1, data2):
	"""
	calculate the socres for peak acc/vel/dis.
	score = S(max|data1|, max|data2|) 
	"""
	update()
	p1 = np.amax(np.absolute(data1))
	p2 = np.amax(np.absolute(data2))
	score = S(p1, p2)
	return p1, p2, score 

def I(data, dt):
	# I(t) = max|integral(data^2)dt|
	return np.amax(np.cumsum(data*data)*dt)

def cal_SI(data1, data2, dt):
	"""
	score = S(IA1, IA2) for Arias intensity 
	score = S(IE1, IE2) for Energy integral 
	"""
	update()
	I1 = I(data1, dt)
	I2 = I(data2, dt)
	SI = S(I1, I2)
	return I1, I2, SI 


def N(data, dt):
	"""N = Ie(t)/IE = Ia(t)/IA"""
	return np.cumsum(data*data)*dt/I(data, dt)

def F(N1, N2):
	return np.absolute(N1-N2)


def cal_SD(data1, data2, dt):
	"""SD = 10*(1-max(F))"""
	update()
	N1 = N(data1, dt)
	N2 = N(data2, dt)
	SD = 10*(1-np.amax(F(N1, N2)))

	return SD 


def cal_Sfs(signal1, signal2, fmin, fmax):
	"""
	calculate the score for Fourier Spectra
	Sfs = mean(S(FS1, FS2))
	"""
	update()
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
	update()
	x = np.cumsum(a1*a2)*dt
	y = np.cumsum(a1*a1)*dt
	z = np.cumsum(a2*a2)*dt
	c = x/(np.power(y, 0.5)*np.power(z, 0.5))
	cc = 10*np.amax(c, 0)
	cc = abs(cc)
	return cc 


def cal_Ssa(signal1, signal2, fmin, fmax):
	"""Calculate the score for Response Spectra"""
	update()
	period = get_period(1/fmax, 1/fmin)
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
	update()

	D1 = duration(signal1)
	D2 = duration(signal2)
	return D1, D2, S(D1, D2)



#  ========================================================= GENERATING =======================================================

def scores_matrix(station1, station2, bands):
	"""
	Generate the 3D matrix of scores 
	"""
	print "...Generating main matrix..."
	bands.insert(0, bands[len(bands)-1])
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = c11 = avg = 0.0

	matrix = np.empty((4, len(bands)+1, 13))
	parameter = np.empty((3, 12))

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


			dt = signal1.dt 

			c1 = cal_SD(signal1.accel, signal2.accel, dt)
			c2 = cal_SD(signal1.velo, signal2.velo, dt)

			# parameter1, parameter2, score for intensity
			a1, a2, c3 = cal_SI(signal1.accel, signal2.accel, dt)
			e1, e2, c4 = cal_SI(signal1.velo, signal2.velo, dt)


			# parameter1, parameter2, score for peak data 
			pga1, pga2, c5 = cal_peak(signal1.accel, signal2.accel)
			pgv1, pgv2, c6 = cal_peak(signal1.velo, signal2.velo)
			pgd1, pgd2, c7 = cal_peak(signal1.displ, signal2.displ)

			c8 = cal_Ssa(signal1, signal2, fmin, fmax)
			c9 = cal_Sfs(signal1, signal2, fmin, fmax)
			c10 = cal_C(signal1.accel, signal2.accel, signal1.dt)

			# duration1, duration2, score
			d1, d2, c11 = cal_D(signal1, signal2)

			scores = np.array([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11], float)
			T = ((c1+c2)/2 + (c3+c4)/2 + c5 + c6 + c7 + c8 + c9 + c10 + c11)/9 
			A = (c1+c2+c3+c4+c5+c6+c7+c8+c9+c10)/10 
			scores = np.insert(scores, 0, T)
			scores = np.insert(scores, 1, A)
			scores = np.around(scores, decimals=2)

			matrix[i][j] = scores 

			# getting parameters used to calculate peak, AI, EI, and duration 
			# for broad band only 
			if j == 1: 
				parameter[i-1] = np.array([pgd1, pgd2, pgv1, pgv2, pga1, pga2, a1, a2, e1, e2, d1, d2], float)

		SA = np.array([],float)
		CA = np.array([],float)

		# calculate the average score of all bands 
		for j in range(0, 13):
			# SA = avg(B1...Bn)
			avg2 = np.average(matrix[i][:,j][1:len(bands)-1])
			SA = np.append(SA, avg2)
			SA = np.around(SA, decimals=2)

			# CA = avg(BB...Bn)
			avg1 = np.average(matrix[i][:,j][:len(bands)-1])
			CA = np.append(CA, avg1)
			CA = np.around(CA, decimals=2)

		matrix[i][-1] = CA
		matrix[i][-2] = SA


	# insert the slide contain all AVERAGE values in front 
	for i in range(0, len(bands)+1):
		for j in range(0, 13):
			average = (matrix[1][i][j] + matrix[2][i][j] + matrix[3][i][j])/3
			matrix[0][i][j] = round(average, 2)
	
	return parameter, matrix

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

def parameter_to_list(parameter):
	"""convert the parameter matrix to list"""
	p = []
	for i in range(0, 12):
		para = parameter[:,i]
		# get the maximum of peak values 
		if 0 <= i < 6:
			p.append(max(para[0], para[1], para[2]))

		for j in range(0, 3):
			p.append(para[j])
	
	# p = [PGD1, PGDNS1, PGDEW1, PGDUD1, PGD2, PGDNS2, PGDEW2, PGDUD2, PGV1, PGVNS1, PGVEW1....]
	return p 
# end of parameter_to_list

# =================================================================== PRINTING ===================================================================
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

def print_scores(filenames, coord, path, parameter, matrix):
	""" generate the file containing all the scores of a list of files. """

	try:
		f = open(path, 'a')
	except IOError, e:
		print e

	file1 = filenames[0].split('/')[-1]
	file2 = filenames[1].split('/')[-1]

	scores = [file1, file2]
	for i in range(0, len(matrix)):
		for j in range(0, 13):
			# reading matrix slide by column 
			col = matrix[i][:,j]
			scores.append(col[-1]) #CA
			scores.append(col[-2]) #SA
			for k in range(0, len(matrix[i])-2):
				# append BB...Bn
				scores.append(col[k])
	
	# insert coordinates of station and epi_distance
	scores.insert(2, coord[0])
	scores.insert(3, coord[1])
	scores.insert(4, coord[2])

	# if require to print the parameters used to get scores
	if parameter:
		scores = scores[:5] + parameter + scores[5:]

	d = '{:>12}'*2 + '{:>12.2f}'*(len(scores)-2) + '\n'
	f.write(d.format(*scores))
	f.close()
# end of print_scores

def set_labels(bands):
	# generate labels for scores file 
	o = ['A', 'N', 'E', 'U']
	b = ['CA', 'SA'] 
	s = ['T', 'A', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']
	b_label = "BB" 

	for i in range(1, len(bands)):
		b_label += ',B'+str(i)
	b_label = b_label.split(',')
	b += b_label

	labels = ['#SIGNAL1', 'SIGNAL2', 'X_COOR', 'Y_COOR', 'EPI_DIS']
	for i in range(0, len(o)):
		for j in range(0, len(b)):
			for k in range(0, len(s)):
				labels.append(o[i]+'_'+b[j]+'_'+s[k])
	return labels
# end of set_labels

def update_labels(labels):
	# add labels for the parameters used to calculate scores
	o = ['_', '_NS_', '_EW_', '_UD_']
	p = ['PGD', 'PGV', 'PGA', 'A', 'E', 'DUR']
	d = ['D', 'S']
	p_labels = []

	for i in range(0, len(p)):
		for j in range(0, len(d)):
			for k in range(0, len(o)):
				if i >= 3 and k == 0:
					pass 
				else: 
					p_labels.append(p[i] + o[k] + d[j])
	labels = labels[:5] + p_labels + labels[5:]
	return labels
# end of update_labels

