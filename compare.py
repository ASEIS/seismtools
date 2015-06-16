#!/usr/bin/env python
# ==========================================================================
# The program is to get filename\directory in various ways, and then call 
# comparison on files.
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
				filename = f + '/' + fp
				if not filename in file_list: 
					file_list.append(filename)

		# if is an non-empty file 
		elif os.path.isfile(f) and os.stat(f).st_size != 0:
			# if the file is .her file; append to list. save for comparison. 
			if f.endswith(".her"):
				compare_list1.append(f)

			# if the file is .txt file containing acceleration data, save for comparison 
			elif f.endswith(".txt"):
				try:
					fp = open(f, 'r')
				except IOError, e:
					print e
					pass 
				lines = fp.read().split('\n')
				if '#' and 'Samples:' and 'dt:' in lines[0]: 
					compare_list2.append(f)
				else: 
					print "[ERROR]: unable to recognize file type."

			# if the file contians a list of files
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
		else:
			print "[ERROR]: no such file or directory."
			return 

def compare1(compare_list):
	"""
	The function is to call comparison between .HER files 

	"""
	compare_list = list(set(compare_list)) # sort and remove duplicated filenames in list 
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
	compare_list = list(set(compare_list)) # sort and remove duplicated filenames in list 
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


get_filename()
read_list(file_list) 

compare1(compare_list1)
# fmin fmax 
# 0.05 - 5
# plot from 0- fend 
compare2(compare_list2)


# TODO: update exception 