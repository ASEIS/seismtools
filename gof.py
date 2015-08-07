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
from gof_engine import * 
from gof_data_sim import * 

np.seterr(divide='ignore', invalid='ignore')

def get_epicenter():
	"""get x and y coordinates of epicenter"""
	epi = ''
	x = 0.0
	y = 0.0 
	while not epi:
		epi = raw_input('== Enter the X and Y coordinates of epicenter: ')
		epi = epi.replace(',', ' ')
		epi = epi.split()
		if len(epi) == 2:
			try: 
				x = float(epi[0])
				y = float(epi[1])
				return x, y 
			except ValueError:
				print "[ERROR]: invalid coordinates."
				epi = ''

		print "[ERROR]: invalid coordinates."
		epi = ''
# end of get_epicenter


def get_in():
	""" get the path of input directories """
	indir1 = ''
	indir2 = ''

	while not indir1:
		indir1 = raw_input('== Enter name of 1st input directory: ')

	while not indir2:
		indir2 = raw_input('== Enter name of 2nd input directory: ')

	# check the existence of two directories
	if (not os.path.exists(indir1)) or (not os.path.exists(indir2)):
		print "[ERROR]: input directory does not exist."
		return get_in()

	return indir1, indir2


def get_out():
	""" get the path of output directory and output file from user"""
	outdir = ''
	outname1 = ''
	outname2 = ''

	# get the destination saving outputs 
	while not outdir: 
		outdir = raw_input('== Enter name of the directory to store outputs: ')

	# check existence of target directory 
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	while not outname1:
		outname1 = raw_input('== Enter name of scores file: ')

	while not outname2:
		outname2 = raw_input('== Enter name of matrics file: ')

	path1 = outdir + '/' + outname1
	path2 = outdir + '/' + outname2
	return  path1, path2


def get_files(): 
	file1 = ''
	file2 = ''
	filelist = ''
	# list1 = list2 = np.array([],str)
	list1 = []
	list2 = []
	coorX = []
	coorY = []

	if len(sys.argv) == 2:
		filelist =  sys.argv[1]

	elif len(sys.argv) == 3: 
		file1 = sys.argv[1]
		file2 = sys.argv[2]

	# if receive two files from user 
	if file1 and file2:
		return file1, file2

	# if receive a file containing a list of files 
	if filelist:
		try:
			# list1, list2 = np.genfromtxt(filelist, comments = '#', dtype = 'str', unpack = True)
			# return list1.tolist(), list2.tolist()
			f = open(filelist, 'r')
		except IOError:
			print "[ERROR]: error loading filelist."
			return False

		for line in f: 
			if not ('#' in line): 
				line = line.split()
				# not containing coordinate 
				if len(line) == 2: 
					list1.append(line[0])
					list2.append(line[1])
					coorX.append(0.0)
					coorY.append(0.0)

				# containing coordinate
				elif len(line) == 4: 
					list1.append(line[0])
					list2.append(line[1])
					try: 
						coorX.append(float(line[2]))
						coorY.append(float(line[3]))
					except ValueError:
						coorX.append(0.0)
						coorY.append(0.0)

		return list1, list2, coorX, coorY


	# if encounter other inputs
	print "[ERROR]: Invalid inputs."
	return False
# end of get_files

def search_file(dirname, info):
	"""search for files contains given network code and station name"""
	file_dir = {'HN':'', 'V1':'', 'BH':'', 'V2':''}
	for fp in os.listdir(dirname):
		if (info in fp) and ('HN' in fp): 
			file_dir['HN'] = fp
		elif (info in fp) and ('V1' in fp): 
			file_dir['V1'] = fp
		elif (info in fp) and ('BH' in fp): 
			file_dir['BH'] = fp 
		elif (info in fp) and ('V2' in fp): 
			file_dir['V2'] = fp
	return file_dir
# end of search_file

def get_bands():
	"""
	The function is to allow user specify sample rates. 
	Without user input, sample rates are setting to default values. 
	"""
	f0 = 0.05
	f1 = 0.1
	f2 = 0.25 
	f3 = 0.5
	f4 = 1
	f5 = 2
	f6 = 4 
	bands = [f0, f1, f2, f3, f4, f5, f6]
	freq = []
	flag = True 

	while flag: 
		flag = False 
		freq = raw_input('== Enter the sequence of sample rates: ').replace(',', ' ').split()
		if not freq:
			#setting to default values 
			return bands 

		if len(freq) == 1:
			print "[ERROR]: invalid sample rates"
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
	return bands 
# enf of get_bands

# ================================================================== READING ================================================================
def read_stamp(filename):
	"""get the time stamp from file's header"""
	try: 
		with open(filename) as f: 
			try: 
				header = f.readlines()[0].split()
				stamp = header[4].split(',')[-1].split(':')
				tmp = stamp[2].split('.')
				stamp[2] = tmp[0]
				stamp.append(tmp[1])

				f.close()
			except IndexError:
				print "[ERROR]: missing time stamp."
				return []
	except IOError:
		print "[ERROR]: No such file."
		return []

	# converting time stamps to floats
	for i in range(0, len(stamp)):
		stamp[i] = float(stamp[i])
	return stamp 
# end of read_stamp


def read_file(filename):
	"""
	The function is to read 10-column .her files. 
	Return a list of psignals for each orientation. 
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
	return station 
# end of read_file

# =========================================================== MAIN =================================================================
def check_data(station):
	"""checks the data after rotation, process_dt, and synchronization 
	to avoid encountering errors in gof_engine """
	for i in range(0, len(station)):
		signal = station[i]

		if signal.accel.size == 0: 
			print "[ERROR]: Empty array after processing signals."
			return False 
		if signal.velo.size == 0: 
			print "[ERROR]: Empty array after processing signals."
			return False 
		if signal.displ.size == 0: 
			print "[ERROR]: Empty array after processing signals."
			return False 

		if np.isnan(np.sum(signal.accel)):
			print "[ERROR]: NaN data after processing signals."
			return False 
		if np.isnan(np.sum(signal.velo)):
			print "[ERROR]: NaN data after processing signals."
			return False 
		if np.isnan(np.sum(signal.displ)):
			print "[ERROR]: NaN data after processing signals."
			return False 
	return station
# end of check_data 

def main(file1, file2):
	"""reads two files and gets two stations; 
	then processes signals in each stations 
	station1 = data 
	station2 = simulation """

	station1 = read_file(file1)
	station2 = read_file(file2)

	if (not station1) or (not station2):
		return False, False 

	station2 = rotate(station2) # rotate simulation 

	# process signals to have the same dt 
	station1, station2 = process_dt(station1, station2)

	# print station1[0].samples
	# print station1[0].accel.size 
	# station1[0].print_attr()
	# station2[0].print_attr()

	# synchronize starting and ending time of data arrays 
	stamp = read_stamp(file1) # get time stamp from data file 
	station1, station2 = synchronize(station1, station2, stamp)


	# station1[0].print_attr()
	# station2[0].print_attr()

	if station1[0].samples != station2[0].samples: 
		print "[ERROR]: two files do not have the same number of samples after processing."
		return False, False 
		# print station1[0].accel.size
		# print station2[0].accel.size

	station1 = check_data(station1)
	station2 = check_data(station2)


	return station1, station2
# end of main

if __name__ == "__main__":
	# getting files or lists of files from user; return tuple 
	files = get_files()

	if isinstance(files[0], str):
		file1, file2 = files
		coor = []

		path = get_out()
		bands = get_bands()

		station1, station2 = main(file1, file2)
		if station1 and station2:
			parameter, matrix = scores_matrix(station1, station2, coor, bands)
			print_matrix(path, matrix)
		else: 
			pass 

	elif isinstance(files[0], list):
		list1, list2, coorX, coorY = files
		indir1, indir2 = get_in()

		# if coordinates are given; ask for epicenter 
		if coorX[0] and coorY[0]:
			Ex, Ey = get_epicenter()
		else: 
			Ex = Ey = 0.0 

		# s_path = get_out('scores')
		# m_path = get_out('matrics')
		s_path, m_path = get_out()

		bands = get_bands()

		try:
			f = open(s_path, 'w')
			m = open(m_path, 'w')
		except IOError, e:
			print e

		labels = set_labels(bands)
		m_labels = set_mlabels()

		d = '{:>12}'*2 + '{:>12.8}'*(len(labels)-2) + '\n'
		# d = '{:>12}'*len(labels) + '\n'
		f.write(d.format(*labels))

		d = '{:>12}'*2 + '{:>12.6}'*(len(m_labels)-2) + '\n'
		m.write(d.format(*m_labels))
		f.close()
		m.close()

		i = 0 
		# for i in range(0, len(list1)):
		while i < len(list1): 
			file1 = indir1 + '/' + list1[i]
			file2 = indir2 + '/' + list2[i]

			# if file1 does not exist, search for files 
			if not os.path.isfile(file1):
				file_dir = search_file(indir1, list1[i])
				for fp in file_dir.values():
					# if found matched file; update the lists 
					if fp: 
						list1.insert(i+1, fp)
						list2.insert(i+1, list2[i])
						coorX.insert(i+1, coorX[i])
						coorY.insert(i+1, coorY[i])

			else: 
				print "[...processing " + file1 + ' and ' + file2 + '...]'

				x = coorX[i]
				y = coorY[i]

				epdist = math.sqrt((x-Ex)**2+(y-Ey)**2)
				coord = [x, y, epdist]

				station1, station2 = main(file1, file2)
				if station1 and station2: 
					parameter, matrix = scores_matrix(station1, station2, bands)
					parameter = parameter_to_list(parameter)

					print_scores([file1,file2], coord, s_path, [], matrix) # print scores 
					print_scores([file1,file2], coord, m_path, parameter, np.array([])) # print values used to calculate scores
				else: 
					pass 
			i += 1 

	print "[DONE]"
	
