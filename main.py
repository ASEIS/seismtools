#!/usr/bin/env python
# ==========================================================================
# The program is to get filename\directory in various ways, and then use 
# the files to call other programs.
# ==========================================================================
import sys
import os
from smc import *
from compare_signals import *

file_list = []
compare_list1 = []
compare_list2 = []

def get_filename(): 
	global file_list
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
	global compare_list1
	global compare_list2

	for f in file_list:
		# if is a directory; read all files\directories in the directory and append them to file_list
		if os.path.isdir(f):
			for fp in os.listdir(f):
				file_list.append(f + '/' + fp)

		# if is a file 
		elif os.path.isfile(f):
			# print f 
			# if the file is V1/raw data file: generate a list of records, text file for acceleration, and .her file 
			if f.endswith(".V1") or f.endswith(".raw"):
				record_list, network, station_id = load_smc_v1(f)
				filename = network + "." + station_id + ".V1.her"
				precord_list = process_record_list(network, station_id, record_list)
				if precord_list == False:
					print_unprocessed(f)
					# append the filename to unprocessed files 
				else: 
					print_her(filename, precord_list) 

			# if the file is V2/processed data file; generate a list of precords, text file for acceleration, and .her file 
			elif f.endswith(".V2"):
				record_list, network, station_id = load_smc_v2(f)
				filename = network + "." + station_id + ".V2.her"
				if print_her(filename, record_list) == False:
					print_unprocessed(f)
					# append the filename to unprocessed files 

			# if the file is .her file; append to list. save for comparison. 
			elif f.endswith(".her"):
				compare_list1.append(f)
			else: 
				# fp = open()
				try:
					fp = open(f, 'r')
				except IOError, e:
					print e
					pass 
				ftype = ' '
				for line in fp:
					# if the file contains a list of filenames 
					if '#filelist' in line: ftype = 'list'
					# if the file contains acceleration data 
					elif '#' and 'Samples:' and 'dt:' in line: 
						ftype = 'txt' 
						break 
					elif ftype == 'list': 
						file_list = file_list + line.split()
					else: 
						print "[ERROR]: unable to recognize file type."
						break 
				if ftype == 'txt': # append to list. save for comparison 
					compare_list2.append(f)



		else:
			print "[ERROR]: no such file or directory."
			return 
	# return a list of her files and a list of txt files that may be used for comparison 
	# return compare_list1, compare_list2

def compare1(compare_list):
	"""
	The function is to call comparison between .HER files 

	"""
	while compare_list:
		for f1 in compare_list:
			for f2 in compare_list:
				# EXAMPLE: CI.Q0028.V1.her v.s. CI.Q0028.V2.her 
				if f2 != f1 and f2[0:-5] == f1[0:-5]: 
					compare_her(f1, f2)
					# compare_list.remove(f1)
					compare_list.remove(f2)
					break 
			compare_list.remove(f1)
			# print compare_list

def compare2(compare_list):
	"""
	The function is to call comparison between .TXT files 
	""" 
	while compare_list:
		for f1 in compare_list:
			for f2 in compare_list:
				# EXAMPLE: CI.Q0028.V1N.txt v.s. CI.Q0028.V2N.txt 
				if f2 != f1 and f2[:-6] == f1[:-6] and f1[-5:] == f2[-5:]: 
					compare_txt(f1, f2)
					# compare_list.remove(f1)
					compare_list.remove(f2)
					break 
			compare_list.remove(f1)
			# print compare_list

def print_unprocessed(filename):
    """
    The function generates a file containing a list of files that were not processed by this program. 
    """
    # generate a text file (header + data)
    f = open('unprocessed_files.txt', 'a')
    f.write(filename +"\n")
    f.close()
# end of print_unprocessed 

# get_filename()
# read_list(file_list) 
# print file_list

# compare1(compare_list1)
# compare2(compare_list2)


# TODO: update exception 