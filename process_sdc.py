#!/usr/bin/env python
# ===================================================================================
# The program is to get filename\directory in various ways, and then processes
# found ASCII files from same station with same sample rate --> generate .her files. 
# ===================================================================================
import sys
import os
from sdc import *

destination = ''
file_list = []
path = ''
event = ''
net = ''
station = ''
info = ''

def split(filename):
	"""
	Split the filename to get path/event.net.station.info.ASCII
	"""
	global path 
	global event 
	global net 
	global station 
	global info 
	
	tmp = filename.split('/')[-1]
	path = filename.replace('/'+tmp, '')
	tmp = tmp.split('.')

	# tmp = [event, net, stat, info]
	if len(tmp) >= 1:
		event = tmp[0]
	if len(tmp) >= 2: 
		net = tmp[1]
	if len(tmp) >= 3:
		station = tmp[2]
	if len(tmp) >= 4:
		info = tmp[3]
		if len(info) == 3:
			info = info[:-1]

# end of split  

def get_station():
	"""To get the event ID; network code; and station ID from user as a searching option"""
	global event 
	global net 
	global station 

	event = raw_input('== Enter event ID (optional): ')
	net = raw_input('== Enter network code (optional): ').upper()
	station = raw_input('== Enter station ID (optional): ').upper()

	# return event, net, station 


def get_info():
	"""
	To get the information code (representing sample rate and data type)
	from user as searching options.
	"""
	global info 
	info = raw_input('== Enter sample rate and data type representation (optional): ').upper()
	if len(info) == 3:
		info = info[:-1]

	# return info 

def search(): 
	"""
	The function is to get a list of files to search for their pairs. 
	"""
	global destination
	global path 
	global event 
	global net 
	global station 
	global info 
	path_list = []


	# if filenam is not given with command 
	if len(sys.argv) == 1:
		while not path_list: 
			path_list = raw_input('== Enter the file\directory path: ')
		path_list = path_list.split()

	# one or more filenames are given; get a list of them
	else: 
		path_list = sys.argv[1:]


	# print path_list
	# iterate through paths; search for files with matched options 
	for p in path_list:
		if not os.path.isdir(p):
			split(p)

		else:
			path = p 

		if not path:
			print "[ERROR]: missing directory path."
			return 

		if (not event) or (not net) or not station:
			get_station()

		if not info:
			get_info()

		for fp in os.listdir(path):
			if (event in fp) and (net in fp) and (station in fp) and ('.'+info in fp) and (fp.endswith('.ascii')):
				file_list.append(fp)

	

	while not destination: 
		destination = raw_input('== Enter name of the directory to store outputs: ')
	# check existence of target directory 
	if not os.path.exists(destination):
		os.makedirs(destination)

	get_destination(destination)

	return file_list

# end of search 



def search_pairs(file_list):
	"""
	search for the pairs in a given list of files. 
	"""
	if not file_list:
		return 

	file_list = sorted(file_list)
	orientation = ['N', 'E', 'Z']
	for f in file_list:
		file_dict = {}
		tmp = f.split('/')[-1]
		tmp = tmp.split('.')
		if len(tmp) < 5:
			return 
		info = tmp[3]
		info = info[0:2]
		for i in range(0, 3):
			# target = tmp[0] + '.' + tmp[1] + '.' + tmp[2] + '.' + info + orientation[i] + '.ascii'
			target = '.' + tmp[1] + '.' + tmp[2] + '.' + info + orientation[i]
			for ff in file_list:
				if target in ff:
					file_dict[orientation[i]] = path + '/' + ff
					if not ff == f: 
						file_list.remove(ff) 

			# if target in file_list:
			# 	file_dict[orientation[i]] = path + '/' + target
			# 	file_list.remove(target)


		# print file_dict

		# process the directory with sdc.py
		if len(file_dict) == 3: 
			print_her(file_dict)
		else:
			print "[ERROR]: No pair found."
		
	# print file_list 
# end of search_pairs
# ------------------------------------------------------------------------------------------------------------------------------

file_list = search()
search_pairs(file_list)





