#!/usr/bin/env python
# ==========================================================================
# The program is to get filename\directory in various ways, and then processes
# found V1/V2 files --> generate column acceleration .txt files and .her files. 
# ==========================================================================
import sys
import os
from smc import *

destination = ''

def get_filename(): 
	file_list = []
	global destination

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
	else: 
		# clear content of unprocessed and warning files if they exist 
		if os.path.exists(destination + '/unprocessed.txt'):
			clear(destination + '/unprocessed.txt')
		if os.path.exists(destination + '/warning.txt'):
			clear(destination + '/warning.txt')

	get_destination(destination)
	return file_list
# end of get_filename


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

				# if encounters errors with records in station 
				if not station.list: 
					station = False 
				else: 
				# process records in station 
					station.process()

				if not station: 
					print_message(f, 'unprocessed')
				else: 
					print_smc(station)
					print_her(station)
					check_station(station)

			# if the file is V2/processed data file; generate text file for acceleration, and .her file 
			elif f.upper().endswith(".V2"):
				station = load_smc_v2(f)

				if not station: 
					print_message(f, 'unprocessed')
				else:
					print_smc(station)
					print_her(station)
					check_station(station)

			else:
				try:
					fp = open(f)
				except IOError, e:
					print e
					pass
				lines = fp.read().split()
				if not '#filelist' in lines[0]:
					print "[ERROR]: unable to recognize file type: " + f
				else:
					for l in lines[1:]:
						if not l in file_list:
							file_list.append(l)
						# file_list += lines[1:]
		else:
			print "[ERROR]: no such file or directory: " + f
# end of read_list 


def print_message(message, ftype):
	"""
	The function is to generate a files containing warning/unprocessed messages for input files. 
	"""
	f = open(destination + '/' + ftype + '.txt', 'a')
	f.write(message +"\n")
	f.close()
# end of print_message


def check_station(station):
    """
    The function is to check the station name of each record,
    if it's in the location should be discarded, print warning. 
    """
    # check instance 
    if not isinstance(station, seism_station):
        return 
    if not station.list:
        return 

    discard = {'dam': 'Dam', 'Fire Sta': 'Fire Station', 'Acosta Res': 'Acosta Res', 'Bldg': 'Building', 'Br': 'Interchange Bridge'}
    name = station.name 
    for key in discard:
        if key in name:
        	filename = station.network + station.id + '.' + station.type 
        	msg = filename + " was processed, but it's from " + discard[key]
        	print_message(msg, 'warning')
        	break 
# end of check_station


def clear(filename):
	# clear the content of a file if it exists 
	try: 
		open(filename, 'w').close()
	except IOError, e:
		pass
# end of clear 



file_list = get_filename()
read_list(file_list) 


