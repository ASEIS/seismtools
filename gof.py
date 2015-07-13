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
import time

np.seterr(divide='ignore', invalid='ignore')

f0 = 0.05
f1 = 0.1
f2 = 0.25 
f3 = 0.5
f4 = 1
f5 = 2
f6 = 4 
bands = [f0, f1, f2, f3, f4, f5, f6]

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
	outname = ''

	# get the destination saving outputs 
	while not outdir: 
		outdir = raw_input('== Enter name of the directory to store outputs: ')

	# check existence of target directory 
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	while not outname:
		outname = raw_input('== Enter name of score file: ')

	path = outdir + '/' + outname
	return  path 


def get_files(): 
	file1 = ''
	file2 = ''
	filelist = ''
	list1 = list2 = np.array([],str)

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
			# list1, list2 = np.loadtxt(filelist, comments = '#', unpack = True)
			list1, list2 = np.genfromtxt(filelist, comments = '#', dtype = 'str', unpack = True)
			return list1.tolist(), list2.tolist()
		except IOError:
			print "[ERROR]: error loading filelist."
			return False, False

	# if encounter other inputs
	print "[ERROR]: Invalid inputs."
	return False, False
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

# ================================================================== READING ================================================================
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

if __name__ == "__main__":
	# getting files or lists of files from user; return tuple 
	files = get_files()

	toolbar_width = 40

	# setup toolbar
	sys.stdout.write("[%s]" % (" " * toolbar_width))
	sys.stdout.flush()
	sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

	for i in xrange(toolbar_width):
	    time.sleep(0.1) # do real work here
	    # update the bar
	    sys.stdout.write("-")
	    sys.stdout.flush()

	sys.stdout.write("\n")

	if isinstance(files[0], str):
		file1, file2 = files
		station1 = read_file(file1)
		station2 = read_file(file2)

		path = get_out()
		get_bands()

		matrix = scores_matrix(station1, station2)
		print_matrix(path, matrix)

	elif isinstance(files[0], list):
		list1, list2 = files
		indir1, indir2 = get_in()
		path = get_out()
		get_bands()

		try:
			f = open(path, 'w')
		except IOError, e:
			print e

		labels = set_labels()
		d = '{:>12}'*len(labels) + '\n'
		f.write(d.format(*labels))
		f.close()

		for i in range(0, len(list1)):
			file1 = indir1 + '/' + list1[i]
			file2 = indir1 + '/' + list2[i]
			station1 = read_file(file1)
			station2 = read_file(file2)
			matrix = scores_matrix(station1, station2)
			print_scores(path, matrix)
	


	# print_matrix(path, matrix)


	# else:
	# 	return 


# for i in range(0, 4):
# 	for j in range(0, len(bands)+1):
# 		print matrix[i][j]
# 	print "---------------------------------------------------------------------------------------"


