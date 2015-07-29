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


def load_event(eventfile):
	"""The function is to read the event file and get information about the event."""
	if not eventfile.lower().endswith('evnt'):
		print "[ERROR]: load event file only. "
		return 
	try:
		fp = open(eventfile)
	except IOError, e:
		print e
		pass 
	header = fp.read().split()
	print header 
	event_id = header[0]
	tmp = header[2].split(',')
	date = tmp[0]
	time = tmp[1]

	latitude = header[3]
	longitude = header[4]
# end of load_event





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
			try: 
				dt = float(tmp[6])
			except ValueError:
				pass 
	
			break 


	if instr_type in data_type:
		dtype = data_type[instr_type]

	signal = seism_signal(samples, dt, data, dtype)
	header = "# " + network + " " + station + " " + "ASCII" + " " + date + " " + str(samples) + " " + str(dt) + "\n"

	return signal
# end of load_file


def process(signal):
	"""
	The function takes a signal, use its's current data to get acceleration, velocity, and displacement. 
	Then return a psignal. 
	"""
	if not isinstance(signal, seism_signal):
		print "[ERROR]: instance error; process signal objects only. "
		return 

	acc = np.array([],float)
	vel = np.array([],float)
	dis = np.array([],float)
	dt = signal.dt 

	if signal.type == 'a':
		
		acc = signal.data 
		window = taper('all', 20, signal.samples)
		acc = acc*window

		vel = integrate(acc, dt)
		vel = s_filter(vel, signal.dt, type = 'highpass', family = 'ellip')

		dis = integrate(vel, dt)
		dis = s_filter(dis, signal.dt, type = 'highpass', family = 'ellip')

	elif signal.type == 'v':
		vel = signal.data 
		window = taper('all', 20, signal.samples)
		vel = vel*window


		acc = derivative(vel, dt)
		acc = s_filter(acc, signal.dt, type = 'highpass', family = 'ellip')

		dis = integrate(vel, dt)
		dis = s_filter(dis, signal.dt, type = 'highpass', family = 'ellip')


	elif signal.type == 'd':
		dis = signal.data 
		window = taper('all', 20, signal.samples)
		dis = dis*window

		vel = derivative(dis, dt)
		vel = s_filter(vel, signal.dt, type = 'highpass', family = 'ellip')

		acc = derivative(vel, dt)
		acc = s_filter(acc, signal.dt, type = 'highpass', family = 'ellip')

	else:
		pass 

	psignal = seism_psignal(signal.samples, signal.dt, np.c_[dis, vel, acc], 'c', acc, vel, dis)
	return psignal
# end of process


def print_her(file_dict):
    """
    The function generates .her files for each station (with all three channels included)
    """
    global destination
     # if there are more than three channels, save for later 
    if len(file_dict) > 3:
        print "==[The function is processing files with 3 channels only.]=="
        return False 

    # compose filename
    filename = file_dict['N'].split('/')[-1]
    filename = filename.replace('N.ascii', '.her')

    try:
        f = open(destination + '/' + filename, 'w')
    except IOError, e:
        print e

    # load files in dictionary; generate siganls and processes them
    file_dict['N'] = process(load_file(file_dict['N']))
    file_dict['E'] = process(load_file(file_dict['E']))
    file_dict['Z'] = process(load_file(file_dict['Z']))

    dis_ns = []
    vel_ns = []
    acc_ns = []
    dis_ew = []
    vel_ew = []
    acc_ew = []
    dis_up = []
    vel_up = []
    acc_up = []

    for key in file_dict:
        if key == 'N':
            dis_ns = file_dict[key].displ.tolist()
            vel_ns = file_dict[key].velo.tolist()
            acc_ns = file_dict[key].accel.tolist()
        elif key == 'E':
            dis_ew = file_dict[key].displ.tolist()
            vel_ew = file_dict[key].velo.tolist()
            acc_ew = file_dict[key].accel.tolist()
        elif key == 'Z':
            dis_up = file_dict[key].displ.tolist()
            vel_up = file_dict[key].velo.tolist()
            acc_up = file_dict[key].accel.tolist()

    signal = file_dict[key]
    # get a list of time incremented by dt 
    time = [0.000]
    samples = file_dict['N'].samples 
    dt = file_dict['N'].dt
    
    while samples > 1:
        time.append(time[len(time)-1] + dt)
        samples -= 1 
    
    network = filename.split('.')[1]
    station = filename.split('.')[2]
    info = filename.split('.')[3]


    # TODO: get time and date from read_event()
    # header = "# " + network + " " + station + " " + "ASCII" + " " + "date" + "," + "time" + " " + str(signal.samples) + " " + str(signal.dt) + "\n"
    f.write(header)

    descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    f.write(descriptor.format("# time", "dis_ns", "dis_ew", "dis_up", "vel_ns", "vel_ew", "vel_up", "acc_ns", "acc_ew", "acc_up")) # header 

    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up):
        f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
    f.close()
    print "*Generated .her file at: " + destination + "/" + filename
#end of print_her 
