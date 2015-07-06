#!/usr/bin/env python
# ================================================================================
# The program is to call comparison on TWO 10-column files or TWO 1-column files. 
# ================================================================================
import sys
import os
# from smc import *
from compare_signals import *
# from stools import *

def check_type(filename):
	"""checks the type of file being .her or .txt"""
	try: 
		with open(filename) as f:
			content = f.readlines()
			num_col = len(content[1].split())
			if num_col == 1:
				return 'TXT'
			elif num_col == 10:
				return 'HER'
			else:
				print "[ERROR]: Cannot recognize file type."
				return False 

	except IOError, e:
		print e 
		return False  
# end of check_type

def get_filename():
	file1 = ''
	file2 = ''

	# get paths of two files 
	if len(sys.argv) >= 2:
		file1 = sys.argv[1]

	if len(sys.argv) >= 3:
		file2 = sys.argv[2]

	while not file1:
		file1 = raw_input('== Enter the path of file1: ')

	while not file2:
		file2 = raw_input('== Enter the path of file2: ')

	if check_type(file1) == check_type(file2) == 'TXT':
		compare_txt(file1, file2)
	elif check_type(file1) == check_type(file2) == 'HER':
		compare_her(file1, file2)
	else:
		return 

# end of get_filename

def get_bound():
	fmin = raw_input('== Enter fmin for FAS: ')
	fmax = raw_input('== Enter fmax for FAS: ')
	try: 
		fmin = float(fmin)
		fmax = float(fmax)
	except ValueError:
		print "[ERROR]: invalid input type: floats or integers ONLY."
		return get_bound()

	while fmin >= fmax:
		print "[ERROR]: fmax must be greater than fmin."
		return get_bound()

	# set_bound(fmin, fmax)
	return fmin, fmax 

get_filename()
