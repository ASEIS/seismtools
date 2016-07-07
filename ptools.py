#!/usr/bin/env python
"""
# ===================================================================================
# The program contains general functions what may be used by other programs.
# ===================================================================================
"""
from __future__ import division, print_function, absolute_import
import os
import sys
import numpy as np
from seism import s_filter, seism_psignal
from stools import seism_cutting, seism_appendzeros

def synchronize_all_stations(obs_data, stations, stamp, eqtimestamp, leading):
    """
    synchronize the stating time and ending time of data arrays
    obs_data = recorded data (optional); stations = simulation signal(s)
    """
    # If we have a recorded data time stamp
    if stamp is not None and obs_data is not None:
        start = stamp[0]*3600 + stamp[1]*60 + stamp[2]
        eq_time = eqtimestamp[0]*3600 + eqtimestamp[1]*60 + eqtimestamp[2]
        sim_start = eq_time - leading

        for i in range(0, 3):
            # synchronize the start time
            if start < sim_start:
                # data time < sim time < earthquake time; cutting data array
                obs_data[i] = seism_cutting('front', (sim_start - start),
                                            20, obs_data[i], False)
            elif start > eq_time:
                # sim time < earthquake time < data time; adding zeros in front
                obs_data[i] = seism_appendzeros('front', (start - eq_time),
                                                20, obs_data[i])
                for station in stations:
                    station[i] = seism_cutting('front', (eq_time - sim_start),
                                               20, station[i], False)
            else:
                # sim time < data time < earthquake time; adding zeros
                obs_data[i] = seism_appendzeros('front', (start - sim_start),
                                                20, obs_data[i])

    # synchronize the ending time
    if obs_data is not None:
        obs_dt = obs_data[0].dt
        obs_samples = obs_data[0].samples
        obs_time = obs_dt * obs_samples
    else:
        obs_time = None

    # Find target timeseries duration
    target_time = None
    if obs_time is not None:
        target_time = obs_time
    for station in stations:
        station_dt = station[0].dt
        station_samples = station[0].samples
        station_time = station_dt * station_samples
        if target_time is None:
            target_time = station_time
            continue
        target_time = min(target_time, station_time)

    # Work on obs_data
    if obs_data is not None:
        for i in range(0, 3):
            if obs_time > target_time:
                obs_data[i] = seism_cutting('end', (obs_time - target_time),
                                            20, obs_data[i], False)
        obs_samples = obs_data[0].samples
        obs_time = obs_dt * obs_samples

    # Work on simulated data
    for station in stations:
        for i in range(0, 3):
            sim_dt = station[i].dt
            sim_samples = station[i].samples
            sim_time = sim_dt * sim_samples
            if sim_time > target_time:
                station[i] = seism_cutting('end', (sim_time - target_time),
                                           20, station[i], False)

    # scale the data if they have one sample in difference after synchronizing
    total_samples = None
    if obs_data is not None:
        total_samples = obs_samples
    for station in stations:
        sim_samples = station[0].samples
        if total_samples is None:
            total_samples = sim_samples
            continue
        total_samples = max(sim_samples, total_samples)

    # For obs_data
    if obs_data is not None:
        for i in range(0, 3):
            if obs_data[i].samples == total_samples - 1:
                obs_data[i] = seism_appendzeros('end', obs_data[i].dt,
                                                20, obs_data[i])
    # For simulated data
    for station in stations:
        for i in range(0, 3):
            if station[i].samples == total_samples - 1:
                station[i] = seism_appendzeros('end', station[i].dt,
                                               20, station[i])

    return obs_data, stations
# end of synchronize_all_stations

def filter_data(psignal, fmin, fmax):
    """
    This function is used to filter with a bandpass filter between fmin/fmax
    """
    if not isinstance(psignal, seism_psignal):
        print("[ERROR]: found error filtering psignal.")
        return False
    delta_t = psignal.dt
    psignal.accel = s_filter(psignal.accel, delta_t, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)
    psignal.velo = s_filter(psignal.velo, delta_t, type='bandpass',
                            family='butter', fmin=fmin, fmax=fmax,
                            N=4, rp=0.1, rs=100)
    psignal.displ = s_filter(psignal.displ, delta_t, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)

    psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

    return psignal
# end of filter_data

def get_files():
    """
    This function parses the parameters specified by the user in the
    command-line. If there's a single argument we treat it as a file
    containing a list of files to be processed. If there are two
    arguments we treat them as two input files to be compared.
    """
    file1 = ''
    file2 = ''
    filelist = ''
    list1 = []
    list2 = []
    coor_x = []
    coor_y = []

    if len(sys.argv) == 2:
        filelist = sys.argv[1]
    elif len(sys.argv) == 3:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    else:
        print("Error: please provide two files to compare or file list!")
        sys.exit(-1)

    # if received two files from user
    if file1 and file2:
        return file1, file2

    # if received a file containing a list of files
    if filelist:
        try:
            input_file = open(filelist, 'r')
        except IOError:
            print("[ERROR]: error loading filelist.")
            return False

        for line in input_file:
            if not '#' in line:
                line = line.split()

                if len(line) == 2:
                    # not containing coordinates
                    list1.append(line[0])
                    list2.append(line[1])
                    coor_x.append(0.0)
                    coor_y.append(0.0)

                elif len(line) == 4:
                    # containing coordinates
                    list1.append(line[0])
                    list2.append(line[1])
                    try:
                        coor_x.append(float(line[2]))
                        coor_y.append(float(line[3]))
                    except ValueError:
                        coor_x.append(0.0)
                        coor_y.append(0.0)

        # Close the input file
        input_file.close()

        return list1, list2, coor_x, coor_y

    # if encounter other inputs
    print("[ERROR]: Invalid inputs.")
    return False
# end of get_files

def get_bands():
    """
    The function is to allow user specify sample rates.
    Without user input, sample rates are setting to default values.
    """
    freq_0 = 0.05
    freq_1 = 0.1
    freq_2 = 0.25
    freq_3 = 0.5
    freq_4 = 1
    freq_5 = 2
    freq_6 = 4
    bands = [freq_0, freq_1, freq_2, freq_3, freq_4, freq_5, freq_6]
    freqs = []
    flag = True

    while flag:
        flag = False
        freqs = raw_input('== Enter the sequence of '
                          'sample rates: ').replace(',', ' ').split()
        if not freqs:
            #setting to default values
            return bands

        if len(freqs) == 1:
            print("[ERROR]: invalid sample rates")
            flag = True
        else:
            bands = []
            for freq in freqs:
                try:
                    bands.append(float(freq))
                except ValueError:
                    print("[ERROR]: invalid sample rates")
                    flag = True
                    break

            for i in range(0, len(bands)-1):
                if bands[i] >= bands[i+1]:
                    print("[ERROR]: invalid sequence of sample rates")
                    flag = True
                    break
    return bands
# end of get_bands

def get_output_format():
    """
    This function asks the user to select between the bbp or her output formats
    """
    output_format = ''

    # Get the output format the user wants
    while output_format != 'bbp' and output_format != 'her':
        output_format = raw_input('== Enter output format (bbp/her): ')
        output_format = output_format.lower()

    return output_format
#end of get_output_format

def check_data(station):
    """
    Checks the data after rotation, process_dt, and synchronization
    to avoid encountering errors in gof_engine
    """
    for i in range(0, len(station)):
        signal = station[i]

        if signal.accel.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if signal.velo.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if signal.displ.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if np.isnan(np.sum(signal.accel)):
            print("[ERROR]: NaN data after processing signals.")
            return False
        if np.isnan(np.sum(signal.velo)):
            print("[ERROR]: NaN data after processing signals.")
            return False
        if np.isnan(np.sum(signal.displ)):
            print("[ERROR]: NaN data after processing signals.")
            return False
    return station
# end of check_data

# ================================ READING ================================
def read_file(filename):
    """
    This function reads a timeseries file(s) in either bbp
    format (ends with .bbp) or hercules (ends otherwise)
    """
    if filename.lower().endswith(".bbp"):
        # Filename in bbp format
        return read_file_bbp(filename)
    # Otherwise use hercules format
    return read_file_her(filename)
# end of read_file

def read_file_bbp2(filename):
    """
    This function reads a bbp file and returns the timeseries in the
    format time, n/s, e/w, u/d tuple
    """
    time = np.array([])
    ns_comp = np.array([])
    ew_comp = np.array([])
    ud_comp = np.array([])

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            line = line.strip()
            if line.startswith('#') or line.startswith('%'):
                # Skip comments
                continue
            # Trim in-line comments
            if line.find('#') > 0:
                line = line[:line.find('#')]
            if line.find('%') > 0:
                line = line[:line.find('%')]
            # Make them float
            pieces = line.split()
            pieces = [float(piece) for piece in pieces]
            time = np.append(time, pieces[0])
            ns_comp = np.append(ns_comp, pieces[1])
            ew_comp = np.append(ew_comp, pieces[2])
            ud_comp = np.append(ud_comp, pieces[3])
    except IOError:
        print("[ERROR]: error reading bbp file.")
        return np.array([]), np.array([]), np.array([]), np.array([])

    # All done!
    return time, ns_comp, ew_comp, ud_comp
# end of read_file_bbp2

def read_file_bbp(filename):
    """
    This function reads timeseries data from a set of BBP files
    """
    # Get filenames for displacement, velocity and acceleration bbp files
    work_dir = os.path.dirname(filename)
    base_file = os.path.basename(filename)

    base_tokens = base_file.split('.')[0:-2]
    dis_tokens = list(base_tokens)
    vel_tokens = list(base_tokens)
    acc_tokens = list(base_tokens)

    dis_tokens.append('dis')
    vel_tokens.append('vel')
    acc_tokens.append('acc')

    dis_tokens.append('bbp')
    vel_tokens.append('bbp')
    acc_tokens.append('bbp')

    dis_file = os.path.join(work_dir, '.'.join(dis_tokens))
    vel_file = os.path.join(work_dir, '.'.join(vel_tokens))
    acc_file = os.path.join(work_dir, '.'.join(acc_tokens))

    # Read 3 bbp files
    [time, dis_ns, dis_ew, dis_up] = read_file_bbp2(dis_file)
    [_, vel_ns, vel_ew, vel_up] = read_file_bbp2(vel_file)
    [_, acc_ns, acc_ew, acc_up] = read_file_bbp2(acc_file)

    samples = dis_ns.size
    delta_t = time[1]

    # samples, dt, data, acceleration, velocity, displacement
    psignal_ns = seism_psignal(samples, delta_t, np.c_[dis_ns, vel_ns, acc_ns],
                               'c', acc_ns, vel_ns, dis_ns)
    psignal_ew = seism_psignal(samples, delta_t, np.c_[dis_ew, vel_ew, acc_ew],
                               'c', acc_ew, vel_ew, dis_ew)
    psignal_up = seism_psignal(samples, delta_t, np.c_[dis_up, vel_up, acc_up],
                               'c', acc_up, vel_up, dis_up)

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of read_file_bbp

def read_file_her(filename):
    """
    The function is to read 10-column .her files.
    Return a list of psignals for each orientation.
    """
    time, dis_ns, dis_ew, dis_up = [np.array([], float) for _ in xrange(4)]
    vel_ns, vel_ew, vel_up = [np.array([], float) for _ in xrange(3)]
    acc_ns, acc_ew, acc_up = [np.array([], float) for _ in xrange(3)]

    try:
        (time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew,
         vel_up, acc_ns, acc_ew, acc_up) = np.loadtxt(filename,
                                                      comments='#',
                                                      unpack=True)
    except IOError:
        print("[ERROR]: error loading her file.")
        return False

    samples = dis_ns.size
    delta_t = time[1]

    # samples, dt, data, acceleration, velocity, displacement
    psignal_ns = seism_psignal(samples, delta_t, np.c_[dis_ns, vel_ns, acc_ns],
                               'c', acc_ns, vel_ns, dis_ns)
    psignal_ew = seism_psignal(samples, delta_t, np.c_[dis_ew, vel_ew, acc_ew],
                               'c', acc_ew, vel_ew, dis_ew)
    psignal_up = seism_psignal(samples, delta_t, np.c_[dis_up, vel_up, acc_up],
                               'c', acc_up, vel_up, dis_up)

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of read_file_her

def read_unit_bbp(filename):
    """
    Get the units from the file's header
    Returns either "m" or "cm"
    """
    units = None

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("units=") > 0:
                units = line.split()[2]
                break
        input_file.close()
    except IOError:
        print("[ERROR]: No such file.")
        sys.exit(-1)

    # Make sure we got something
    if units is None:
        print("[ERROR]: Cannot find units in bbp file!")
        sys.exit(-1)

    # Figure out if we have meters or centimeters
    if units == "cm" or units == "cm/s" or units == "cm/s^2":
        return "cm"
    elif units == "m" or units == "m/s" or units == "m/s^2":
        return "m"

    # Invalid units in this file
    print("[ERROR]: Cannot parse units in bbp file!")
    sys.exit(-1)
# end of read_unit_bbp

def read_stamp(filename):
    """
    Get the time stamp from file's header
    """
    if filename.endswith(".bbp"):
        # File in bbp format
        return read_stamp_bbp(filename)
    # Otherwise use hercules format
    return read_stamp_her(filename)
# end of read_stamp

def read_stamp_bbp(filename):
    """
    Get the time stamp from the bbp file's header
    """
    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("time=") > 0:
                stamp = line.split()[2].split(',')[-1].split(':')
                break
        input_file.close()
    except IOError:
        print("[ERROR]: No such file.")
        return []

    # Converting time stamps to floats
    stamp = [float(i) for i in stamp]
    return stamp
# end of read_stamp_bbp

def read_stamp_her(filename):
    """
    Get the time stamp from the her file's header
    """
    try:
        with open(filename) as input_file:
            try:
                header = input_file.readline().split()
                stamp = header[4].split(',')[-1].split(':')
                input_file.close()
            except IndexError:
                print("[ERROR]: missing time stamp.")
                return []
    except IOError:
        print("[ERROR]: No such file.")
        return []

    # converting time stamps to floats
    for i in range(0, len(stamp)):
        stamp[i] = float(stamp[i])
    return stamp
# end of read_stamp_her

# ================================ WRITING ==================================
def print_her(filename, station):
    # filename = 'processed-' + filename.split('/')[-1]
    try:
        out_f = open(filename, 'w')
    except IOError, e:
        print(e)
    dis_ns = station[0].displ.tolist()
    vel_ns = station[0].velo.tolist()
    acc_ns = station[0].accel.tolist()
    dis_ew = station[1].displ.tolist()
    vel_ew = station[1].velo.tolist()
    acc_ew = station[1].accel.tolist()
    dis_up = station[2].displ.tolist()
    vel_up = station[2].velo.tolist()
    acc_up = station[2].accel.tolist()

    # get a list of time incremented by dt
    time = [0.000]
    samples = station[0].samples
    dt = station[0].dt
    tmp = samples

    while tmp > 1:
        time.append(time[len(time)-1] + dt)
        tmp -= 1

    out_f.write('# missing header \n')

    descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    out_f.write(descriptor.format("# time",
                                  "dis_ns", "dis_ew", "dis_up",
                                  "vel_ns", "vel_ew", "vel_up",
                                  "acc_ns", "acc_ew", "acc_up")) # header

    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time,
                                                      dis_ns, dis_ew, dis_up,
                                                      vel_ns, vel_ew, vel_up,
                                                      acc_ns, acc_ew, acc_up):
        out_f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9))
    out_f.close()
# end of print_her

def print_bbp(input_file, output_file, station):
    """
    This function generates processed .bbp files for
    each of velocity/acceleration/displacement
    and copies the header of the input bbp file
    """
    output_dir = os.path.dirname(output_file)
    output_basename = os.path.basename(output_file)

    # Prepare data for output
    acc_ns = station[0].accel.tolist()
    vel_ns = station[0].velo.tolist()
    dis_ns = station[0].displ.tolist()
    acc_ew = station[1].accel.tolist()
    vel_ew = station[1].velo.tolist()
    dis_ew = station[1].displ.tolist()
    acc_up = station[2].accel.tolist()
    vel_up = station[2].velo.tolist()
    dis_up = station[2].displ.tolist()

    # Start with time = 0.0
    time = [0.000]
    samples = station[0].samples
    while samples > 1:
        time.append(time[len(time)-1] + station[0].dt)
        samples -= 1

    # Prepare to output
    out_data = [['dis', dis_ns, dis_ew, dis_up, 'displacement', 'cm'],
                ['vel', vel_ns, vel_ew, vel_up, 'velocity', 'cm/s'],
                ['acc', acc_ns, acc_ew, acc_up, 'acceleration', 'cm/s^2']]

    for data in out_data:
        if not output_basename.endswith('.bbp'):
            # Remove extension
            bbp_output_basename = os.path.splitext(output_basename)[0]
            bbp_output_filename = os.path.join(output_dir,
                                               "%s.%s.bbp" %
                                               (bbp_output_basename,
                                                data[0]))
            output_header = ["# Station: NoName",
                             "#    time= 00/00/00,00:00:00.00 UTC",
                             "#     lon= 0.00",
                             "#     lat= 0.00",
                             "#   units= %s" % (data[5]),
                             "#",
                             "# Data fields are TAB-separated",
                             "# Column 1: Time (s)",
                             "# Column 2: N/S component ground "
                             "%s (+ is 000)" % (data[4]),
                             "# Column 3: E/W component ground "
                             "%s (+ is 090)" % (data[4]),
                             "# Column 4: U/D component ground "
                             "%s (+ is upward)" % (data[4]),
                             "#"]
        else:
            # Read header of input file
            input_dirname = os.path.dirname(input_file)
            input_basename = os.path.basename(input_file)
            pieces = input_basename.split('.')
            pieces = pieces[0:-2]
            bbp_input_file = os.path.join(input_dirname,
                                          "%s.%s.bbp" %
                                          ('.'.join(pieces),
                                           data[0]))
            input_header = []
            in_fp = open(bbp_input_file, 'r')
            for line in in_fp:
                line = line.strip()
                if line.startswith("#"):
                    input_header.append(line)
            in_fp.close()

            # Compose new header
            output_header = []
            for item in input_header:
                if item.find("units=") > 0:
                    output_header.append("#   units= %s" % (data[5]))
                else:
                    output_header.append(item)

            pieces = output_basename.split('.')
            pieces = pieces[0:-2]
            bbp_output_filename = os.path.join(output_dir,
                                               "%s.%s.bbp" %
                                               ('.'.join(pieces),
                                                data[0]))
        # Write output file
        try:
            out_fp = open(bbp_output_filename, 'w')
        except IOError, e:
            print(e)
            continue

        # Write header
        for item in output_header:
            out_fp.write("%s\n" % (item))

        # Write timeseries
        for val_time, val_ns, val_ew, val_ud in zip(time, data[1],
                                                    data[2], data[3]):
            out_fp.write("%5.7f   %5.9e   %5.9e    %5.9e\n" %
                         (val_time, val_ns, val_ew, val_ud))

        # All done, close file
        out_fp.close()
        print("*Writing file: %s " % (bbp_output_filename))
# end of print_bbp
