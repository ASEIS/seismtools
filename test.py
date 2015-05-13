#!/usr/bin/env python
import os.path
filename = "NC.NHC.E.txt"
num = 0 
while os.path.isfile(filename):
	num += 1
	filename = filename[:-3] + str(num) + filename[-4:]
	print filename
	num += 1
	filename = filename[:-3] + str(num) + filename[-4:]
	print filename
	num += 1
	filename = filename[:-3] + str(num) + filename[-4:]
	print filename

# NC.NHC.E.1.txt
# NC.NHC.E.1.2.txt
# NC.NHC.E.1.2.3.txt
