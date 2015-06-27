#!/usr/bin/env python
# ===================================================================================
# The program is to get filename\directory in various ways, and then processes
# found ASCII files from same station with same sample rate --> generate .her files. 
# ===================================================================================
import sys
import os
# from sdc import *

destination = ''
file_list = []
event = ''
net = ''
station = ''

def get_option():
	"""
	The function is to get evernt ID; station; and information code (representing sample rate and data type)
	from user as searching options.
	"""
	event = raw_input('== Enter event ID (optional): ')
	network = raw_input('== Enter network code (optional): ').upper()
	station = raw_input('== Enter station ID (optional): ').upper()
	info = raw_input('== Enter sample rate and data type representation (optional): ').upper()

	return event, network+'.'+station, info 


	pass 
# end of get_option

def get_path(): 
	"""
	The function is to get a list of files to search for their pairs. 
	"""
	global search_list
	global destination

	# if filenam is not given with command 
	if len(sys.argv) == 1:
		while not search_list: 
			path_list = raw_input('== Enter the file\directory path: ')
		path_list = path_list.split()

	# one or more filenames are given; get a list of them
	else: 
		path_list = sys.argv[1:]

	# iterate through paths; search for files with options 
	for path in path_list:
		if os.path.isdir(path):
			event, station, info = get_option()
			filename = path + event + '.' + station + info 
			search(filename)

		else:
			search(path)

	while not destination: 
		destination = raw_input('== Enter name of the directory to store outputs: ')
	# check existence of target directory 
	if not os.path.exists(destination):
		os.makedirs(destination)
	get_destination(destination)

# end of get_path

def search_dir(path):
	"""Return all files in directory."""
	file_list = []
	for fp in os.listdir(path):
		file_list.append(fp)
	return sorted(file_list)


def search_station(file_list):
	"""Return a list of data files from given station."""
	global event 
	global net 
	global station 
	target = event + '.' + net + '.' + station
	l = []

	for f in file_list:
		if target in f: 
			l.append(f)
			file_list.remove(f)
	return sorted(l)


def search_info(file_list):
	global info 
	if len(info) == 3:
		info = info.replace(info[2], '')

	for i in range(0, 3):
		info += orientation[i]
		# f = event + '.' + net + '.' + station + '.' + info  
		for f in file_list:
			if info in f:
				l.append(f)
				file_list.remove(f)
		info = info.replace(info[2], '')
	return sorted(l)






	# 	for i in range(0, 3):
	# 		info += orientation[i]
	# 		f = event + '.' + net + '.' + station + '.' + info 
	# 		print f 
	# 		for fp in os.listdir(path):
	# 			if f in fp: 
	# 				print fp 
	# 				file_list.append(fp)
	# 		info = info.replace(info[2], '')
	# 	return file_list

	# pass 




def search(filename):
	"""
	The function is to search for the pairs of given files. 
	"""

	file_dict = {}
	file_list = []
	orientation = ['N', 'E', 'Z']

	tmp = filename.split('/')[-1]
	path = filename.replace(tmp, '')
	tmp = tmp.split('.')

	# tmp = [event, net, stat, info]
	if len(tmp) < 4:
		print "[ERROR]: invalid filename/path."
		return 

	event = tmp[0]
	net = tmp[1]
	station = tmp[2]

	# if containing information code for searching 
	if tmp[3]:
		info = tmp[3]
		if len(info) == 3:
			info = info.replace(info[2], '')


		for i in range(0, 3):
			info += orientation[i]
			f = event + '.' + net + '.' + station + '.' + info 
			print f 
			for fp in os.listdir(path):
				if f in fp: 
					print fp 
					file_list.append(fp)
			info = info.replace(info[2], '')
		return file_list

	
	# if containing network and station id for searching 
	elif tmp[1] and tmp[2]:
		f = event + '.' + net + '.' + station 
		for fp in os.listdir(path):
			if f in fp: 
				file_list.append(fp)
		file_list = sorted(file_list)
		print file_list 
		return file_list

	# i = len(tmp) - 2
	# for i in range(0, len(tmp)):
	# 	filename = tmp[0] + 
	# tmp 
	
	# tmp = tmp.split('.')
	# if 'ascii' in tmp: tmp.remove('ascii')
	# if len(tmp) >= 4: 
	# 	info = tmp[3]
	# 	if len(info) == 3: 
	# 		file_list.append()

	# event = tmp[0]
	# network = tmp[1].upper()
	# station_id = tmp[2].upper()
	# info = tmp[3]

	# for s in search_list:
	# 	# if is a directory
	# 	if os.path.isdir(s):
	# 		for fp in os.listdir(s):
	# 			filename = s + '/' + fp
	# 			if not filename in file_list: 
	# 				file_list.append(filename)

	# 	# if is an existing file 
	# 	elif os.path.isfile(s) and os.stat(f).st_size != 0:



	# 	tmp = s.split('/')[-1]
	# 	path = s.replace(tmp, "")
	# 	print path 
	# 	tmp = tmp.split('.')
	# 	event = tmp[0]
	# 	network = tmp[1].upper()
	# 	station_id = tmp[2].upper()
	# 	info = tmp[3]
	


	# 	print tmp 




	pass 
# end of search 





# get_search_list()
search('data-sdc/14383980.CI.CHN.')
# read_list(file_list) 


