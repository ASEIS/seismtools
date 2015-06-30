#!/usr/bin/env python
# ===================================================================================
# The program is to read two .her files; process their signals; 
# calculate their scores with different sample rates; 
# and generate 3D matrix for scores. 
# ===================================================================================

import sys
import os
import numpy as np
from seism import *


destination = 'outputs'
bb = 0.05
b1 = 0.1
b2 = 0.25 
b3 = 0.5
b4 = 1
b5 = 2
b6 = 4 
rate = [bb, b1, b2, b3, b4, b5, b6]


def get_files(): 
	file1 = ''
	file2 = ''
	global destination

	# get paths of two files 
	if len(sys.argv) == 2:
		file1 =  sys.argv[1]

	elif len(sys.argv) == 3: 
		file1 = sys.argv[1]
		file2 = sys.argv[2]

	while not file1:
		file1 = raw_input('== Enter the path of file1: ')

	while not file2:
		file2 = raw_input('== Enter the path of file2: ')

	get_rate()


	# get the destination saving outputs 
	# while not destination: 
	# 	destination = raw_input('== Enter name of the directory to store outputs: ')

	# # check existence of target directory 
	# if not os.path.exists(destination):
	# 	os.makedirs(destination)

	# get_destination(destination)



	return file1, file2
# end of get_files

def get_rate():
	"""
	The function is to allow user specify sample rates. 
	Without user input, sample rates are setting to default values. 
	"""
	global rate
	input_rate = []
	flag = True 

	while flag: 
		flag = False 
		input_rate = raw_input('== Enter the sequence of 6 sample rates: ').split()
		if not input_rate:
			#setting to default values 
			return

		if len(input_rate) != 7: 
			print "[ERROR]: missing sample rates."
			flag = True 
		else: 

			try: 
				bb = float(input_rate[0])
				b1 = float(input_rate[1])
				b2 = float(input_rate[2])
				b3 = float(input_rate[3])
				b4 = float(input_rate[4])
				b5 = float(input_rate[5])
				b6 = float(input_rate[6])
			except ValueError: 
				print "[ERROR]: invalid sample rates"
				flag = True 

			if not (bb  < b1 < b2 < b3 < b4 < b5 < b6):
				print "[ERROR]: invalid sequence of sample rates"
				flag = True 

	rate = [bb, b1, b2, b3, b4, b5, b6]
# enf of get_rate

def read_file(filename):
	"""
	The function is to read 10-column .her files. 
	Return a list of psignals for each orientation. 
	"""
	time = dis_ns = dis_ew = dis_up = vel_ns = vel_ew = vel_up = acc_ns = acc_ew = acc_up = np.array([],float)

	try:
		time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, skiprows = 2, unpack = True)
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


def filt(psignal):
	"""The function is to call ellip_filter on each psignal from seism.py"""
	if not isinstance(psignal, seism_psignal):
		print "[ERROR]: encounter error filting psignal."
		return False 

	# order of unlabeled arguments is: (data, dt, order N, rp, rs, Wn) 
	psignal.accel = ellip_filter(psignal.accel, psignal.dt) 
	psignal.velo = ellip_filter(psignal.velo, psignal.dt) 
	psignal.displ = ellip_filter(psignal.displ, psignal.dt) 

	psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

	# psignal.print_attr()

	return psignal 


def S(p1, p2):
	# S(p1, p2) = 10*exp{-[(p1-p2)/min(p1, p2)]^2}
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



def scores_matrix(station1, station2):
	"""
	Generate the 3D matrix of scores 
	"""
	global rate 
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = c11 = avg = 0.0

	matrix = np.empty((3, 8, 12))

	for i in range(0, len(station1)):
		for j in range(0, len(rate)-1): 
			fmin = rate[j]
			fmax = rate[j+1]
			# print fmin, fmax 

			signal1 = filt(station1[i])
			signal2 = filt(station2[i])

			c5 = cal_peak(signal1.accel, signal2.accel)
			c6 = cal_peak(signal1.velo, signal2.velo)
			c7 = cal_peak(signal1.displ, signal2.displ)

			scores = np.array([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11], float)
			scores = np.append(scores, np.average(scores))
			scores = np.around(scores, decimals=2)

			matrix[i][j] = scores 

		fs2 = np.array([],float)
		fs1 = np.array([],float)

		# calculate the average score of all bands 
		# FS2 = avg(b1-b6)
		# FS1 = avg(bb-b6)
		for j in range(0, 12):
			avg2 = np.average(matrix[i][:,j][1:6])
			fs2 = np.append(fs2, avg2)
			fs2 = np.around(fs2, decimals=2)

			avg1 = np.average(matrix[i][:,j][:6])
			fs1 = np.append(fs1, avg1)
			fs1 = np.around(fs1, decimals=2)


		matrix[i][-1] = fs1
		matrix[i][-2] = fs2



		# fs2 = np.append([], avg2)
		# avg1 = np.average(matrix[i][:,0])
		# fs1 = np.append([], avg1)

		# matrix = np.append(matrix, [fs2, fs1])

	
	return matrix

def my_scores_matrix(station1, station2):
	pass 


def print_matrix(matrix):
	"""
	The function generates a text files printing the scores matrix. 
	"""
	global destination
	filename = "scores.txt"
	ll = []
	try:
		f = open(destination + '/' + filename, 'w')
	except IOError, e:
		print e

	for i in range(0, 3): 
		for j in range(0, 8):
			for k in range(0, 12): 
				f.write(str(matrix[i][j][k]) + '\t')
			f.write('\n')

	# 		ll.append(matrix[i][j].tolist())

	# descriptor = '{:>12.2f}'*len(ll) + '\n'

	# 	for i in range(0, len(l)):
	# 		f.write(descriptor.format(l[i]))

		# f.write(descriptor.format(l))
	f.close()
	pass 






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

# print_matrix(matrix)

for i in range(0, 3):
	for j in range(0, 8):
		print matrix[i][j]