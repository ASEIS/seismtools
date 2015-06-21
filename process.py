#!/usr/bin/env python
# ==========================================================================
# The program is to get filename\directory in various ways, and then processes
# found V1/V2 files --> generate column acceleration .txt files and .her files. 
# ==========================================================================
import sys
import os
from smc import *

destination = ''
file_list = []

def get_filename(): 
	global file_list
	global destination
	# clear('unprocessed_files.txt')
	# clear('warning.txt')

	# if filenam is not given with command 
	if len(sys.argv) == 1:
		while not file_list: 
			file_list = raw_input('== Enter the file\directory name: ')
		file_list = file_list.split()

	# one or more filenames are given; get a list of them
	else: 
		file_list = sys.argv[1:]

	while not destination: 
		destination = raw_input('== Enter name of the directory to store outputs: ')
	# check existence of target directory 
	if not os.path.exists(destination):
		os.makedirs(destination)

	get_destination(destination)





def read_list(file_list):
	"""
	The function is to read a list of files/directory and check their 
	types to call corresponding functions. 
	"""

	for f in file_list:
		# if is a directory; read all files\directories in the directory and append them to file_list
		if os.path.isdir(f):
			for fp in os.listdir(f):
				filename = f + '/' + fp
				if not filename in file_list: 
					file_list.append(filename)

		# if is an non-empty file 
		elif os.path.isfile(f) and os.stat(f).st_size != 0:
			# if the file is V1/raw data file: generate text file for acceleration, and .her file 
			if f.upper().endswith(".V1") or f.upper().endswith(".RAW"):
				station = load_smc_v1(f)
				if not station: 
					print_unprocessed(f)
				else: 
					print_smc(station)
					print_her(station)

				# record_list, network, station_id = load_smc_v1(f)
				# filename = network + "." + station_id + ".V1.her"
				# precord_list = process_record_list(network, station_id, record_list)
				# if precord_list == False:
				# 	print_unprocessed(f)
				# 	# append the filename to unprocessed files 
				# else: 
				# 	print_her(filename, precord_list) 

			# if the file is V2/processed data file; generate text file for acceleration, and .her file 
			elif f.upper().endswith(".V2"):
				station = load_smc_v2(f)
				if not station: 
					print_unprocessed(f)
				else:
					print_smc(station)
					print_her(station)
				# record_list, network, station_id = load_smc_v2(f)
				# filename = network + "." + station_id + ".V2.her"
				# if print_her(filename, record_list) == False:
				# 	print_unprocessed(f)
					# append the filename to unprocessed files 

			else: 
				try:
					fp = open(f)
				except IOError, e:
					print e
					pass 
				lines = fp.read().split()
				# print lines 
				if not '#filelist' in lines[0]:
					print "[ERROR]: unable to recognize file type."
				else: 
					for l in lines[1:]: 
						if not l in file_list:
							file_list.append(l)
					# file_list += lines[1:]
		else:
			print "[ERROR]: no such file or directory: " + f
		# file_list = list(set(file_list)) #remove duplicated ones from list 
# end of read_list 

def print_unprocessed(filename):
    """
    The function generates a file containing a list of files that were not processed by this program. 
    """
    f = open(destination + '/unprocessed_files.txt', 'a')
    f.write(filename +"\n")
    f.close()
# end of print_unprocessed 

def clear(filename):
	# clear the content of a file if it exists 
	try: 
		open(filename, 'w').close()
	except IOError, e:
		pass
# end of clear 



get_filename()
read_list(file_list) 


# TODO: update exception 