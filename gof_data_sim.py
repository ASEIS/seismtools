#!/usr/bin/env python
# ===================================================================================
# The program processes data reference and simulation signals. 
# Including rotation, making equal dt, and synchronization of starting and ending time. 
# ===================================================================================
import os
import numpy as np
from scipy import interpolate
from seism import *
from stools import * 
# import matplotlib.pyplot as plt


def get_azimuth():
	""" get the azimuth for rotation from user. """
	azimuth = ''
	while not azimuth: 
		a = raw_input('== Enter azimuth for rotation (optional): ')

		# if user choose not to rotate 
		if not a: 
			return azimuth

		try :
			azimuth = float(a)
		except ValueError:
			print "[ERROR]: invalid azimuth."
	return azimuth

# def rotate(signal1, signal2):
# 	""" Rotation of simulation signal; 
# 	signal1 = data reference
# 	signal2 = simulation signal"""
# 	azimuth = get_azimuth()
# 	if not azimuth:
# 		return 
# 	pass 

def rotate(station):
	""" station = [psignal_ns, psignal_ew, psignal_up] """
	# checking instance 
	if len(station) != 3:
		return 
	for s in station:
		if not isinstance(s, seism_psignal):
			return 
	azimuth = get_azimuth()
	if not azimuth:
		return

	psignal_ns = station[0]
	psignal_ew = station[1]
	psignal_up = station[2]

	# rotate data in North and East 
	matrix = np.array([(math.cos(math.radians(azimuth)), -math.sin(math.radians(azimuth))), (math.sin(math.radians(azimuth)), math.cos(math.radians(azimuth)))])
	[psignal_ns.accel, psignal_ew.accel] = matrix.dot([psignal_ns.accel, psignal_ew.accel])
	[psignal_ns.velo, psignal_ew.velo] = matrix.dot([psignal_ns.velo, psignal_ew.velo])
	[psignal_ns.disp, psignal_ew.disp] = matrix.dot([psignal_ns.displ, psignal_ew.displ])

	# rotate data in Up
	psignal_up.accel *= -1 
	psignal_up.velo *= -1 
	psignal_up.displ *= -1 

	station = [psignal_ns, psignal_ew, psignal_up]
	return station 

# end of rotate 
# =============================================================================================================== 

def get_dt():
	dt = ''
	while not dt: 
		d = raw_input("== Enter common dt of two signals: ")
		try:
			dt = float(d)
		except ValueError:
			print "[ERROR]: invalid dt."
	return dt 
# end of get_dt

def get_fmax():
	fmax = ''
	while not fmax:
		f =  raw_input("== Enter the maximum frequence for decimation: ")
		try:
			fmax = float(f)
		except ValueError:
			print "[ERROR]: invalid fmax."
	return fmax
# end of get_fmax

def interp(data, t, samples, dt):
	""" call interpolate on given data """
	f = interpolate.interp1d(t, data, bounds_error = False)
	new_t = np.arange(0, samples*dt, dt)
	new_data = f(new_t)

	# using plot to test 
	# plt.plot(t,data,'r',new_t,new_data,'b')
	# plt.show()
	return new_data
# end of interpolate

def process(signal, dt, fmax):
	""" processes signal with common dt and fmax."""
	# call low_pass filter at fmax 
	signal.accel = lowpass_filter(signal.accel, dt, fmax)
	signal.velo = lowpass_filter(signal.velo, dt, fmax)
	signal.displ = lowpass_filter(signal.displ, dt, fmax)

	# interpolate 
	t = np.arange(0, signal.samples*signal.dt, signal.dt)
	signal.accel = interp(signal.accel, t, signal.samples, dt)
	signal.velo = interp(signal.velo, t, signal.samples, dt)
	signal.displ = interp(signal.displ, t, signal.samples, dt)

	return signal
# end of process

def process_dt(station1, station2):
	""" process all signals in two stations to have common dt """
	dt = get_dt()
	fmax = get_fmax()

	# process signals in stations 
	for i in range(0, 3):
		station1[i] = process(station1[i], dt, fmax)
		station2[i] = process(station2[i], dt, fmax)
	
	return station1, station2
# end of process_dt

# =============================================================================================================== 



def read_file(filename):
	"""
	The function is to read 10-column .her files. 
	Return a list of psignals for each orientation. 
	"""
	time = dis_ns = dis_ew = dis_up = vel_ns = vel_ew = vel_up = acc_ns = acc_ew = acc_up = np.array([],float)

	try:
		time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up = np.loadtxt(filename, comments='#', unpack = True)
	except IOError:
		print "[ERROR]: error loading her file. "
		return False  

	samples = dis_ns.size 
	dt = time[1]

	# samples, dt, data, acceleration, velocity, displacement 
	psignal_ns = seism_psignal(samples, dt, np.c_[dis_ns, vel_ns, acc_ns], 'c', acc_ns, vel_ns, dis_ns)
	psignal_ew = seism_psignal(samples, dt, np.c_[dis_ew, vel_ew, acc_ew], 'c', acc_ew, vel_ew, dis_ew)
	psignal_up = seism_psignal(samples, dt, np.c_[dis_up, vel_up, acc_up], 'c', acc_up, vel_up, dis_up)

	station = [psignal_ns, psignal_ew, psignal_up]
	return station 
# end of read_file



station1 = read_file('1/2-CICHN-2.sim')
station2 = read_file('1/2-CICHN-1.dat')
process_dt(station1, station2)
# station[0].print_attr()
# station[1].print_attr()
# station[2].print_attr()


# station = rotate(station)
# station[0].print_attr()
# station[1].print_attr()
# station[2].print_attr()
