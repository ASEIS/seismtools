#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from seism import *
# import numpy as np
# import math
# y = 104 
# x = 14
# x = (x+y-90)/2 
# h1 = 20 
# h2 = 90

# # x = 104 
# # y = 14
# # h1 = 20 
# # h2 = 90


# N = h1*math.cos(math.radians(x)) - h2*math.sin(math.radians(x))
# E = h1*math.sin(math.radians(x)) + h2*math.cos(math.radians(x))
# print N 
# print E


# j = np.array([(math.cos(math.radians(x)), -math.sin(math.radians(x))), (math.sin(math.radians(x)), math.cos(math.radians(x)))])
# # w = np.array([h1, h2])
# print(j.dot([h1, h2]))
# print(j.dot([h2, h1]))

# import os
# file_list = raw_input('== Enter the file\directory name: ')
# print len(file_list)
# print type(file_list)


function [f,fs] = fourierbounded(s,fmin,fmax,dt,points);

abs(np.fft.fft(data))*dt 
# fft_s(:,1) = fft(s(:,1),points);
# afs_s = abs(fft_s)*dt; 


freq = (1/dt)*(0:points-1)/points;
# freq = freq';

deltaf = (1/dt)/points;

ini = fix(fmin/deltaf)+1;
fin = fix(fmax/deltaf);

fs = afs_s(ini:fin,:);
f  = freq(ini:fin,:);

def read_file(filename): 
	"""
	The function is to read general 1-column text files. Return a signal object. 
	"""
	try:
		f = open(filename, 'r')
	except IOError, e:
		print e
		return 
	
	samples = 0 
	dt = 0.0 
	data = []
	for line in f:
		# get header 
		if "#" in line: 
			tmp = line.split()
			samples = int(tmp[4])
			dt = float(tmp[6])
		# get data 
		else:
			# print line.split()[0]
			data.append(float(line))
	data = np.array(data)
	f.close()
	return seism_signal(samples, dt, data, 'a')
	# return samples, dt, data 