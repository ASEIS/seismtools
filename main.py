#!/usr/bin/env python
# ==========================================================================
# The program is to get filename\directory in various ways, and then use 
# the files to call other programs.
# ==========================================================================
import sys
import os
from smc import *
from compare_signals import *

# if filenam is not given with command 
if len(sys.argv) == 1:
	file_list = raw_input('== Enter the file\directory name: ')
	file_list = file_list.split()

# one or more filenames are given; get a list of them
else: 
	file_list = sys.argv[1:]

def read_list(file_list):
	"""
	The function is to read a list of files/directory and check their 
	types/contents to call corresponding functions. 
	"""
	new_file_list = []
	compare_list = []
	for f in file_list:
		# if is a directory; read all files\directories in the directory and append them to file_list
		if os.path.isdir(f):
			for fp in os.listdir(f):
				file_list.append(f + '/' + fp)
			print file_list

		# if is a file 
		elif os.path.isfile(f):
			# if the file is V1/raw data file: generate a list of records, text file for acceleration, and .her file 
			if f.endswith(".V1") or f.endswith(".raw"):
				record_list, network, station_id = load_smc_v1(f)
				filename = network + "." + station_id + ".V1.her"
				print_her(filename, process_record_list(network, station_id, record_list)) 

			# if the file is V2/processed data file; generate a list of precords, text file for acceleration, and .her file 
			elif f.endswith(".V2"):
				record_list, network, station_id = load_smc_v2(f)
				filename = network + "." + station_id + ".V2.her"
				print_her(filename, record_list) 

			# if the file is .her file 
			elif f.endswith(".her"):
				compare_list.append(f)
				# TODO: waiting for another her file to call compare()
				pass 

			fp = open(f, 'r')
			ftype = ' '
			for line in fp:
				# if the file contains a list of filenames 
				if '#filelist' in line: ftype = 'list'
				elif ftype == 'list': 
					file_list = file_list + line.split()
				else:
					# TODO: check for acceleration data text file 
					# TODO: waiting for another her file to call compare()
					pass

			print file_list

		else:
			print "[ERROR]: no such file or directory."
			return 


	return new_file_list


newlist = read_list(file_list)
print newlist
