#!/usr/bin/env python
# ==========================================================================
# The program is to load files downloaded with STP client, then create
# signal objects, and generate .her files.
# ==========================================================================
import numpy as np
from seism import *
from stools import *

destination = ''
header = ''
def get_destination(d):
	"""The function is to get the user input from process.py."""
	global destination
	destination = d
# end of get_destination

def load_file(filename):
	"""
	The function is to read general 1-column text files. Return a signal object.
	"""
	if not filename.lower().endswith("ascii"):
		print "[ERROR]: process Waveform files in ascii format only. "
		return

	global header
	data_type = {'H': 'v', 'L': 'v', 'N': 'a'}
	# v for velocity; a for acceleration, d for displacement

	network = ''
	station = ''
	dt = 0.0
	date = ''
	dtype = ''

	f = filename.split('/')[-1]
	tmp = f.split('.')
	event = tmp[0]


	try:
		f = open(filename, 'r')
	except IOError, e:
		print e
		return

	data = np.loadtxt(filename, comments = '#', unpack = True)
	samples = data.size

	for line in f:
		# get header
		if '#' in line:
			tmp = line.split()
			network = tmp[1]
			station = tmp[2]
			info = tmp[3]

			sample_rate = info[0].upper()
			instr_type = info[1].upper()
			orientation = info[2]

			date = tmp[4]
			time = date.split(',')[-1]
			try:
				dt = float(tmp[6])
			except ValueError:
				pass

			break


	if instr_type in data_type:
		dtype = data_type[instr_type]

	signal = seism_signal(samples, dt, data, dtype)
	header = "# " + network + " " + station + " " + "ASCII" + " " + date + " " + str(samples) + " " + str(dt) + "\n"

	return signal, time
# end of load_file


def process(signal):
	"""
	The function takes a signal, use its's current data to get acceleration, velocity, and displacement.
	Then return a psignal.
	"""
	if not isinstance(signal, seism_signal):
		print "[ERROR]: instance error; process signal objects only. "
		return

	# Correction of base lines was commented out here to avoid problems with synchronization
	# correct_baseline(signal)

	acc = np.array([],float)
	vel = np.array([],float)
	dis = np.array([],float)
	dt = signal.dt

	if signal.type == 'a':

		acc = signal.data
		acc = s_filter(acc, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

		vel = integrate(acc, dt)
		vel = s_filter(vel, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

		dis = integrate(vel, dt)
		dis = s_filter(dis, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

	elif signal.type == 'v':

		vel = signal.data
		vel = s_filter(vel, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

		acc = derivative(vel, dt)

		dis = integrate(vel, dt)
		dis = s_filter(dis, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

	elif signal.type == 'd':

		dis = signal.data
		dis = s_filter(dis, signal.dt, type = 'highpass', family = 'butter', fmin = 0.05, N = 5)

		vel = derivative(dis, dt)
		acc = derivative(vel, dt)

	else:
		pass

	psignal = seism_psignal(signal.samples, signal.dt, np.c_[dis, vel, acc], 'c', acc, vel, dis)
	return psignal
# end of process

def synchronize(stamps, signals):
	"""synchronize signals with given time stamps"""
	end_time = []
	start_time = []
	# convert time stamps to time in seconds
	for i in range(0, len(stamps)):
		time = stamps[i]
		signal = signals[i]

		start = time.split(':')
		start = float(start[0])*3600 + float(start[1])*60 + float(start[2])
		end = start + signal.samples*signal.dt

		start_time.append(start) # update time lists
		end_time.append(end)

	start_last = max(start_time)
	end_first = min(end_time)
	index = start_time.index(start_last)
	new_stamp = stamps[index]

	# cutting signals
	for i in range(0, len(signals)):
		start = start_time[i]
		end = end_time[i]
		if start != start_last:
			t_diff = start_last - start
			signals[i] = seism_cutting('front', t_diff, 20, signals[i], True)

		if end != end_first:
			t_diff = end - end_first
			signals[i] = seism_cutting('end', t_diff, 20, signals[i], True)

	if (signals[0].samples == signals[1].samples == signals[2].samples):
		return new_stamp, signals

	# scale the data if they have one sample in difference after cutting
	for i in range(1, len(signals)):
		if signals[i].samples == signals[i-1].samples - 1:
			signals[i-1].data = signals[i-1].data[1:]
			signals[i-1].samples -= 1
		elif signals[i].samples == signals[i-1].samples + 1:
			signals[i].data = signals[i].data[1:]
			signals[i].samples -= 1
	return new_stamp, signals
# end of synchronize

def print_her(file_dict):
	"""
	The function generates .her files for each station (with all three channels included)
	"""
	global destination
	global header
	 # if there are more than three channels, save for later
	if len(file_dict) > 3:
		print "==[The function is processing files with 3 channels only.]=="
		return False

	# compose filename
	filename = file_dict['N'].split('/')[-1]
	filename = filename.replace((filename.split('.')[0]+'.'), '') #remove event ID
	filename = filename.replace('N.ascii', '.her')

	try:
		f = open(destination + '/' + filename, 'w')
	except IOError, e:
		print e

	# load files in dictionary; generate siganls
	signal_ns, time_ns = load_file(file_dict['N'])
	signal_ew, time_ew = load_file(file_dict['E'])
	signal_up, time_up = load_file(file_dict['Z'])

	# correct baselines before synchronizing times to avoid issues with tappers (if any)
	correct_baseline(signal_ns)
	correct_baseline(signal_ew)
	correct_baseline(signal_up)

	# synchronize signals
	if not (signal_ns.samples == signal_ew.samples == signal_up.samples):
		signals = [signal_ns, signal_ew, signal_up]
		stamps = [time_ns, time_ew, time_ew]
		new_stamp, [signal_ns, signal_ew, signal_up] = synchronize(stamps, signals)

		# update header
		tmp = header.split(' ')
		tmp[-2] = str(signal_ns.samples)
		tmp[-3] = tmp[-3].split(',')[0]+','+new_stamp
		# tmp[-3].split(',')[-1] = new_stamp
		header = ''
		for i in range(0, len(tmp)):
			header += tmp[i] + ' '
		# header += '\n'

	# process signals
	signal_ns = process(signal_ns)
	signal_ew = process(signal_ew)
	signal_up = process(signal_up)


	dis_ns = signal_ns.displ.tolist()
	vel_ns = signal_ns.velo.tolist()
	acc_ns = signal_ns.accel.tolist()
	dis_ew = signal_ew.displ.tolist()
	vel_ew = signal_ew.velo.tolist()
	acc_ew = signal_ew.accel.tolist()
	dis_up = signal_up.displ.tolist()
	vel_up = signal_up.velo.tolist()
	acc_up = signal_up.accel.tolist()

	# print len(dis_ns)
	# print len(vel_ns)
	# print len(acc_ns)
	# print len(dis_ew)
	# print len(vel_ew)
	# print len(acc_ew)
	# print len(dis_up)
	# print len(vel_up)
	# print len(acc_up)


	# get a list of time incremented by dt
	time = [0.000]
	samples = signal_ns.samples
	dt = signal_ns.dt
	tmp = samples

	while tmp > 1:
		time.append(time[len(time)-1] + dt)
		tmp -= 1

	network = filename.split('.')[1]
	station = filename.split('.')[2]
	info = filename.split('.')[3]

	f.write(header)

	descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
	f.write(descriptor.format("# time", "dis_ns", "dis_ew", "dis_up", "vel_ns", "vel_ew", "vel_up", "acc_ns", "acc_ew", "acc_up")) # header

	descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
	for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up):
		f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
	f.close()
	print "*Generated .her file at: " + destination + "/" + filename
#end of print_her
