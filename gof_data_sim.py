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


def rotate(station):
	""" station = [psignal_ns, psignal_ew, psignal_up] """
	# checking instance 
	if len(station) != 3:
		return station 
	for s in station:
		if not isinstance(s, seism_psignal):
			return station 

	psignal_ns = station[0]
	psignal_ew = station[1]
	psignal_up = station[2]

	# rotate data in Up/Down
	psignal_up.accel *= -1 
	psignal_up.velo *= -1 
	psignal_up.displ *= -1 


	azimuth = get_azimuth()
	if not azimuth:
		return station

	# rotate data in North and East 
	matrix = np.array([(math.cos(math.radians(azimuth)), -math.sin(math.radians(azimuth))), (math.sin(math.radians(azimuth)), math.cos(math.radians(azimuth)))])
	[psignal_ns.accel, psignal_ew.accel] = matrix.dot([psignal_ns.accel, psignal_ew.accel])
	[psignal_ns.velo, psignal_ew.velo] = matrix.dot([psignal_ns.velo, psignal_ew.velo])
	[psignal_ns.disp, psignal_ew.disp] = matrix.dot([psignal_ns.displ, psignal_ew.displ])

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
	f = interpolate.interp1d(t, data, 'linear', bounds_error = False)
	new_t = np.arange(0, samples*dt, dt)
	new_data = f(new_t)

	# using plot to test 
	# plt.plot(t,data,'r',new_t,new_data,'b')
	# plt.show()

	# print data 
	# print new_data
	return new_data
# end of interpolate

def process(signal, dt, fmax):
	""" processes signal with common dt and fmax."""
	# call low_pass filter at fmax 
	signal.accel = s_filter(signal.accel, dt, type = 'lowpass', family = 'ellip', max = fmax, N = 3, rp = 0.1, rs = 100) 
	signal.velo = s_filter(signal.velo, dt, type = 'lowpass', family = 'ellip', fmax = fmax, N = 3, rp = 0.1, rs = 100) 
	signal.displ = s_filter(signal.displ, dt, type = 'lowpass', family = 'ellip', fmax = fmax, N = 3, rp = 0.1, rs = 100) 

	# interpolate 
	t = np.arange(0, signal.samples*signal.dt, signal.dt)
	signal.accel = interp(signal.accel, t, signal.samples, dt)
	signal.velo = interp(signal.velo, t, signal.samples, dt)
	signal.displ = interp(signal.displ, t, signal.samples, dt)

	signal.dt = dt 

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
def get_earthq():
	"""get the earthquake start time"""
	time =  raw_input("== Enter the earthquake start time (#:#:#.#): ").replace('.', ':')
	time = time.split(':')
	if len(time) < 4:
		print "[ERROR]: invalid time format."
		return get_earthq()
	else: 
		for i in range(0, len(time)):
			try:
				time[i] = float(time[i])
			except ValueError:
				print "[ERROR]: invalid time format."
				return get_earthq()
				
	# time = [hour, min, sec, frac]
	return time 
# end of get_earthq

def get_leading():
	"""get the simulation leading time"""
	lt = ''
	while not lt:
		t = raw_input("== Enter the simulation leading time (sec): ")
		try: 
			lt = float(t)
		except ValueError:
			print "[ERROR]: invalid leading time."
	return lt 
# end of get_leading

def cut_signal(t_diff, signal):
	if not isinstance(signal, seism_psignal):
		return signal
	num = int(t_diff/signal.dt) 
	signal.samples -= num 
	signal.accel = signal.accel[num:]
	signal.velo = signal.velo[num:]
	signal.displ = signal.displ[num:]
	return signal

def add_signal(t_diff, signal):
	if not isinstance(signal, seism_psignal):
		return signal
	num = int(t_diff/signal.dt) 
	zeros = np.zeros(num)
	signal.samples += num 
	
	signal.accel = np.append(zeros, signal.accel)
	signal.velo = np.append(zeros, signal.velo)
	signal.displ = np.append(zeros, signal.displ)

	return signal


def synchronize(station1, station2, stamp):
	"""synchronize the stating time and ending time of data arrays in two signals
	signal1 = data signal; signal2 = simulation signal """
	if not stamp: 
		stamp = [0, 0, 0, 0]
	print stamp

	eq = get_earthq()
	lt = get_leading()

	# time in sec = hr*3600 + min*60 + sec + frac*0.1
	start = stamp[0]*3600 + stamp[1]*60 + stamp[2] + stamp[3]*0.1 
	eq_time = eq[0]*3600 + eq[1]*60 + eq[2] + eq[3]*0.1 
	sim_start = eq_time - lt 

	for i in range(0, 3):
		signal1 = station1[i]
		signal2 = station2[i]

		dt = signal1.dt # same dt of two signals 
		samples = signal1.samples # original samples 

		# synchronize the start time 
		if start < sim_start: 
			# data time < sim time < earthquake time; cutting data array 
			signal1 = cut_signal((sim_start - start), signal1)

		elif start > eq_time:
			# sim time < earthquake time < data time; adding zeros in front 
			signal1 = add_signal((start - eq_time), signal1)
			signal2 = cut_signal((eq_time - sim_start), signal2)

		else: 
			# sim time < data time < earthquake time; adding zeros 
			signal1 = add_signal((start - sim_start), signal1)


		# synchronize the ending time 
		data_time = dt * samples # total time of data signal 
		end = start + data_time
		sim_time = dt * samples # total simulation time 
		sim_end = sim_start + sim_time 

		if sim_end < end: 
			# adding zeros in simulation signal
			num = int((end - sim_end)/dt)
			zeros = np.zeros(num)
			signal2.samples += num 

			signal2.accel = np.append(signal2.accel, zeros)
			signal2.velo = np.append(signal2.velo, zeros)
			signal2.displ = np.append(signal2.displ, zeros)

		elif end < sim_end:
			# cutting from simulation signal 
			num = int((sim_end - end)/dt) 
			signal2.samples -= num 
			num *= -1 
			signal2.accel = signal2.accel[:num]
			signal2.velo = signal2.velo[:num]
			signal2.displ = signal2.displ[:num]
		else:
			pass 


		station1[i] = signal1
		station2[i] = signal2


	return station1, station2
# end of synchronize 

