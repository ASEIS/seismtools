#!/usr/bin/env python
from __future__ import division
import numpy as np
import math 
import os
import sys 
import time 
# z = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], float)
# y = np.array([z, z, z, z, z, z, z, z], float)
# x = np.array([y, y, y], float)



# print x.shape

# data = range(1,5)
# print data 
# print np.average(data[:2])
# print np.average(data)



# period = np.logspace(1/0.05, 1/4, num=20)
# for p in period:
# 	print p 
# print period

# w = 2*math.pi*period
# print w 

# l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# for f in l:
# 	print l
# 	print f 
# 	for f2 in l:
# 		# print l 
# 		if f2 == f*2:
# 			l.remove(f)
# 			l.remove(f2)
# 	print l 
# a = np.log10(1/4)
# b = np.log10(1/0.05) 

# l = np.linspace(a, b, 20)
# l = np.power(10, l)
# print l 

# num_b = 7-2
# r_label = "BB" 
# for i in range(1, num_b):
# 	r_label += ',B'+str(i)
# r_label += ",SA,CA"
# r_label = r_label.split(',')

# print r_label

# descriptor = '{:>12}' + '  {:>12}'*(num_b+1) + '\n'
# print(descriptor.format(*r_label)) 

# def test(i):
# 	if i == 1:
# 		return 1, 2
# 	else:
# 		return '1', '2'
# 	return 

# def test2(i):
# 	if isinstance(test(i)[0], int):
# 		print "is int"
# 	else:
# 		print "is string"

# test2(1)
# test2(2)


# f0 = 0.05
# f1 = 0.1
# f2 = 0.25 
# f3 = 0.5
# f4 = 1
# f5 = 2
# f6 = 4 
# bands = [f0, f1, f2, f3, f4, f5, f6]
# def set_labels():
# 	global bands 
# 	o = ['A', 'N', 'E', 'U']
# 	b = ['CA', 'SA'] 
# 	s = ['T', 'A', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']
# 	b_label = "BB" 
# 	for i in range(1, len(bands)-1):
# 		b_label += ',B'+str(i)
# 	b_label = b_label.split(',')
# 	b += b_label

# 	labels = ['# SIGNAL1'. 'SIGNAL2']
# 	for i in range(0, len(o)):
# 		for j in range(0, len(b)):
# 			for k in range(0, len(s)):
# 				labels.append(o[i]+'_'+b[j]+'_'+s[k])
# 	return labels

# print set_labels()

# widgets = ['Something: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
#            ' ', ETA(), ' ', FileTransferSpeed()]


# from scipy import interpolate
# import matplotlib.pyplot as plt

# x = np.arange(0, 10)
# y = np.exp(-x/3.0)
# f = interpolate.interp1d(x, y)

# xnew = np.arange(0,9, 0.1)
# ynew = f(xnew)   # use interpolation function returned by `interp1d`
# plt.plot(x, y, 'o', xnew, ynew, '-')
# plt.show()

# toolbar_width = 40
# count = 0

# # setup toolbar
# sys.stdout.write("[%s]" % (" " * toolbar_width))
# sys.stdout.flush()
# sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

# for i in xrange(toolbar_width):
# 	count += 1 
# 	time.sleep(0.1) # do real work here

# 	sys.stdout.write("-")
# 	sys.stdout.write(str(count*100/toolbar_width)+'%') 
# 	sys.stdout.flush()

# sys.stdout.write("\n")

# count = '-' 
# sys.stdout.write(count)
# sys.stdout.flush()
# count = ' '
# sys.stdout.write(count)
# sys.stdout.flush()

# import numpy as np
# import matplotlib.pyplot as plt


# # generate some data
# x = np.arange(0, 10, 0.2)
# y = np.sin(x)

# # plot it
# fig = plt.figure(figsize=(8, 6))
# ax0 = plt.subplot2grid((1, 3), (0, 0), colspan=2)
# ax0.plot(x, y)
# ax1 = plt.subplot2grid((1, 3), (0, 2))
# ax1.plot(y, x)

# plt.show()
