#!/usr/bin/env python
# ==========================================================================
# The program is to load files downloaded with STP client, then create 
# signal objects, and generate .her files.
# ==========================================================================
import numpy as np
from seism import *

destination = '' 
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

	band = {'H': '80-250Hz', 'B': '10-80Hz', 'E': '80-250Hz'}
	data_type = {'H': 'v', 'L': 'v', 'N': 'a'}
	# v for velocity; a for acceleration, d for displacement 

	network = ''
	station = ''
	dt = 0.0
	orientation = ''
	date = ''
	dtype = ''

	f = filename.split('/')[-1]
	tmp = f.split('.')
	event = tmp[0]
	network = tmp[1].upper()
	station_id = tmp[2].upper()
	info = tmp[3]

	sample_rate = info[0].upper()
	instr_type = info[1].upper()
	orientation = info[2]

	try:
		f = open(filename, 'r')
	except IOError, e:
		print e
		return 

	data = np.loadtxt(filename, skiprows = 1, unpack = True)
	samples = data.size -1 

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

	if sample_rate in band:
		pass 
		# TODO

	if instr_type in data_type:
		dtype = data_type[instr_type]

	signal = seism_signal(samples, dt, data, dtype)
	filename = event + '.' + network + '.' + station_id + '.' + sample_rate + instr_type
	# record = seism_record(samples, dt, data, dtype, station, '', '', )
	# samples, dt, data, signal type, station, 
 #        location_lati, location_longi, depth, date, time, orientation
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

	if signal.type == 'a':
		acc = signal.data 
		vel = signal.integrate(acc)
		dis = signal.integrate(velocity)

	elif signal.type == 'v':
		vel = signal.data 
		acc = signal.derivative(vel)
		dis = signal.integrate(vel)

	elif signal.type == 'd':
		dis = signal.data 
		vel = signal.derivative(dis)
		acc = signal.derivative(vel)
	else:
		pass 

	psignal = seism_psignal(signal.samples, signal.dt, np.c_[dis, vel, acc], 'a', acc, vel, dis)
	return psignal
# end of process


def print_her(filename, dict):
    """
    The function generates .her files for each station (with all three channels included)
    """
    global destination
     # if there are more than three channels, save for later 
    if len(dict) > 3:
        print "==[The function is processing files with 3 channels only.]=="
        return False 

    try:
        f = open(destination + '/' + filename, 'w')
    except IOError, e:
        print e
        # return 

    dis_ns = []
    vel_ns = []
    acc_ns = []
    dis_ew = []
    vel_ew = []
    acc_ew = []
    dis_up = []
    vel_up = []
    acc_up = []

    for key in dict:
        if key == 'N':
            dis_ns = dict[key].displ.tolist()
            vel_ns = dict[key].velo.tolist()
            acc_ns = dict[key].accel.tolist()
        elif key == 'E':
            dis_ew = dict[key].displ.tolist()
            vel_ew = dict[key].velo.tolist()
            acc_ew = dict[key].accel.tolist()
        elif key == 'Z':
            dis_up = dict[key].displ.tolist()
            vel_up = dict[key].velo.tolist()
            acc_up = dict[key].accel.tolist()

    signal = dict[key]
    # get a list of time incremented by dt 
    time = [0.000]
    samples = dict['N'].samples 
    dt = dict['N'].dt
    
    while samples > 1:
        time.append(time[len(time)-1] + dt)
        samples -= 1 
    
    # header = "# " + station.network + " " + station.id + " " + station.type + " " + precord.date + "," + precord.time + " " + str(precord.samples) + " " + str(precord.dt) + "\n"
    # f.write(header)

    descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    f.write(descriptor.format("time", "dis_ns", "dis_ew", "dis_up", "vel_ns", "vel_ew", "vel_up", "acc_ns", "acc_ew", "acc_up")) # header 

    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up):
        f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
    f.close()
    print "*Generated .her file at: " + destination + "/" + filename
#end of print_her 


# process(load_file('data-sdc/14383980.CI.WNS.BLN.ascii'))
# load_event('data-sdc/14383980.evnt')
# file_list = ['data-sdc/14383980.CI.CHN.HHN.ascii', 'data-sdc/14383980.CI.CHN.HHE.ascii', 'data-sdc/14383980.CI.CHN.HHZ.ascii']



def main(file_list):
	"""
	The function reads a list of files containing data from three orientations. 
	And call load_file and process on each of them. Then print .her file. 
	"""
	dict = {}
	if len(file_list) < 3: 
		print "[ERROR]: Missing file. The program processes three files in a group. "
		return 
		
	dict['N'] = process(load_file(file_list[0]))
	dict['E'] = process(load_file(file_list[1]))
	dict['Z'] = process(load_file(file_list[2]))
	print_her('14383980.CI.WNS.BL.her', dict)

# end of main