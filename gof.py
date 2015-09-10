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
import matplotlib.pyplot as plt

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

# def search_file(dirname, info):
# 	"""search for files contains given network code and station name"""
# 	file_dir = {'HN':'', 'V1':'', 'BH':'', 'V2':''}
# 	for fp in os.listdir(dirname):
# 		if (info in fp) and ('HN' in fp):
# 			file_dir['HN'] = fp
# 		elif (info in fp) and ('V1' in fp):
# 			file_dir['V1'] = fp
# 		elif (info in fp) and ('BH' in fp):
# 			file_dir['BH'] = fp
# 		elif (info in fp) and ('V2' in fp):
# 			file_dir['V2'] = fp
# 	return file_dir
# # end of search_file

def search_file(dirname, info):
	"""
	Search for files contains given station code and name
	"""
	info = info.replace('.', '')
	for fp in os.listdir(dirname):
		tmp = fp.replace('.', '')
		if (info in tmp) and ('HN' in tmp):
			return fp
			# file_dir['HN'] = fp
		elif (info in tmp) and ('V1' in tmp):
			# file_dir['V1'] = fp
			return fp
		elif (info in tmp) and ('BH' in tmp):
			# file_dir['BH'] = fp
			return fp
		elif (info in tmp) and ('V2' in tmp):
			# file_dir['V2'] = fp
			return fp
		# endif case
	# end for

	# was not found, return same info	
	return info
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

# =========================================================== READING ===========================================================
def read_stamp(filename):
	"""get the time stamp from file's header"""
	try:
		with open(filename) as f:
			try:
				header = f.readlines()[0].split()
				stamp = header[4].split(',')[-1].split(':')
				# tmp = stamp[2].split('.')
				# stamp[2] = tmp[0]
				# stamp.append(tmp[1])

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

def print_her(filename, station):
	filename = 'processed-' + filename.split('/')[-1]
	try:
		f = open(filename, 'w')
	except IOError, e:
		print e
	dis_ns = station[0].displ.tolist()
	vel_ns = station[0].velo.tolist()
	acc_ns = station[0].accel.tolist()
	dis_ew = station[1].displ.tolist()
	vel_ew = station[1].velo.tolist()
	acc_ew = station[1].accel.tolist()
	dis_up = station[2].displ.tolist()
	vel_up = station[2].velo.tolist()
	acc_up = station[2].accel.tolist()


	# get a list of time incremented by dt
	time = [0.000]
	samples = station[0].samples
	dt = station[0].dt
	tmp = samples

	while tmp > 1:
		time.append(time[len(time)-1] + dt)
		tmp -= 1

	f.write('# missing header \n')

	descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
	f.write(descriptor.format("# time", "dis_ns", "dis_ew", "dis_up", "vel_ns", "vel_ew", "vel_up", "acc_ns", "acc_ew", "acc_up")) # header

	descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
	for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up):
		f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
	f.close()
# end of print_her

# ========================================================= MAIN ===============================================================
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


def process(station1, station2, azimuth, commondt, decifmax, eq_time, leading):
	"""
	This method processes the signals in each pair of stations.
	Processing consists on scaling, rotation, decimation, alignment
	and other things to make both signals compatible to apply GOF method.
	station 1: data
	station 2: simulaiton
	"""

	# scale synthetics to cm(/s/s)
	station2 = scale_synthetics(station2)

	# rotate synthetics
	station2 = rotate(station2, azimuth)

	# process signals to have the same dt
	station1, station2 = process_dt(station1, station2, commondt, decifmax)

	# Optional plotting for checking
	# signal1 = station1[0]
	# signal2 = station2[0]
	# t1 = np.arange(0, signal1.samples*signal1.dt, signal1.dt)
	# t2 = np.arange(0, signal2.samples*signal2.dt, signal2.dt)
	# plt.plot(t1,signal1.accel,'r',t2,signal2.accel,'b')
	# plt.show()

	# synchronize starting and ending time of data arrays
	stamp = read_stamp(file1) # get time stamp from data file
	station1, station2 = synchronize(station1, station2, stamp, eq_time, leading)

	# Optional plotting for checking
	# signal1 = station1[0]
	# signal2 = station2[0]
	# t = np.arange(0, signal1.samples*signal1.dt, signal1.dt)
	# plt.plot(t,signal1.accel,'r',t,signal2.accel,'b')
	# plt.show()

	if station1[0].samples != station2[0].samples:
		print "[ERROR]: two files do not have the same number of samples after processing."
		return False, False

	station1 = check_data(station1)
	station2 = check_data(station2)

	return station1, station2

# end of process

if __name__ == "__main__":
	
	# getting files or lists of files from user; return tuple
	files = get_files()

	if isinstance(files[0], str):

		# Start: Two-Files Option

		file1, file2 = files
		coor = []

		# captures input data
		s_path, m_path = get_out()
		bands    = get_bands()
		azimuth  = get_azimuth()
		commondt = get_dt()
		decifmax = get_fmax()
		eq_time  = get_earthq()
		leading  = get_leading()

		# reads signals
		station1 = read_file(file1)
		station2 = read_file(file2)

		# processing signals
		if station1 and station2:
			station1, station2 = process(station1, station2, azimuth, commondt, decifmax, eq_time, leading)
		else:
			print "...Ignoring pair:   " + list1[i] + " - " + list2[i]
			pass

		if station1 and station2:
			parameter, matrix = scores_matrix(station1, station2, bands)
			print_matrix(s_path, matrix)
		else:
			pass

		# End: Two-Files Option

	elif isinstance(files[0], list):

		# Start: List of Files Option

		list1, list2, coorX, coorY = files
		indir1, indir2 = get_in()

		# captures input data
		s_path, m_path = get_out()
		bands    = get_bands()
		decifmax = get_fmax()
		commondt = get_dt()
		azimuth  = get_azimuth()
		eq_time  = get_earthq()
		leading  = get_leading()
		if coorX[0] and coorY[0]:
			Ex, Ey = get_epicenter()
		else:
			Ex = Ey = 0.0

		# open output files
		try:
			f = open(s_path, 'w')
			m = open(m_path, 'w')
		except IOError, e:
			print e

		# prepares formats for output files

		labels = set_labels(bands)
		m_labels = set_mlabels()

		d = '{:>12}'*2 + '{:>12.8}'*(len(labels)-2) + '\n'
		f.write(d.format(*labels))

		d = '{:>12}'*2 + '{:>12.6}'*(len(m_labels)-2) + '\n'
		m.write(d.format(*m_labels))
		f.close()
		m.close()

		# loop of the list of pairs given in the list-file
		for i in range(0, len(list1)):

			# capture full path to files
			file1 = indir1 + '/' + list1[i]
			file2 = indir2 + '/' + list2[i]

			# if file1 not in dir1, search for a match
			if not os.path.isfile(file1):
				fp = search_file(indir1, list1[i])
				if fp == list1[i]:
					# if returns without change, move on
					print "...Ignoring pair:   " + list1[i] + " - " + list2[i] + " (no data)"
					continue
				file1 = indir1 + '/' + fp
			# endif

			# if file2 not in dir2, move on
			if not os.path.isfile(file2):
				print "...Ignoring pair:   " + list1[i] + " - " + list2[i] + " (no synthetic)"
				continue
			# endif

			# Both files are available, attempts to process...
			print "\n...Processing pair: " + file1 + " - " + file2

			# computes epicentral distance
			x = coorX[i]
			y = coorY[i]
			epdist = math.sqrt((x-Ex)**2+(y-Ey)**2)
			coord = [x, y, epdist]

			# reads signals
			station1 = read_file(file1)
			station2 = read_file(file2)

			# processing signals
			if station1 and station2:
				station1, station2 = process(station1, station2, azimuth, commondt, decifmax, eq_time, leading)
			else:
				print "...Ignoring pair:   " + list1[i] + " - " + list2[i]
				continue

			# Optional plotting for checking
			# signal1 = station1[1]
			# signal2 = station2[1]
			# t1 = np.arange(0, signal1.samples*signal1.dt, signal1.dt)
			# t2 = np.arange(0, signal2.samples*signal2.dt, signal2.dt)
			# plt.plot(t1,signal1.accel,'r',t2,signal2.accel,'b')
			# plt.show()

			if station1 and station2:
				# print_her(file1, station1)
				# print_her(file2, station2)
				parameter, matrix, flag = scores_matrix(station1, station2, bands)

				# sanity check to avoid division by zero
				if not flag:
					print "\n...Ignoring pair:   " + list1[i] + " - " + list2[i] + " (div by zero)"
					continue
				# end if: sanity check


				parameter = parameter_to_list(parameter)

				print_scores([file1,file2], coord, s_path, [], matrix) # print scores
				print_scores([file1,file2], coord, m_path, parameter, np.array([])) # print values used to calculate scores
			else:
				pass

	print "[DONE]"

