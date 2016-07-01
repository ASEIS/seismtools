#!/usr/bin/env python
"""
# ===================================================================================
# The program contains general functions what may be used by other programs.
# ===================================================================================
"""
from __future__ import division, print_function, absolute_import
import sys
import numpy as np
from seism import s_filter, seism_psignal

def filter_data(psignal, fmin, fmax):
    """
    This function is used to filter with a bandpass filter between fmin/fmax
    """
    if not isinstance(psignal, seism_psignal):
        print("[ERROR]: encounter error filting psignal.")
        return False
    dt = psignal.dt
    psignal.accel = s_filter(psignal.accel, dt, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)
    psignal.velo = s_filter(psignal.velo, dt, type='bandpass',
                            family='butter', fmin=fmin, fmax=fmax,
                            N=4, rp=0.1, rs=100)
    psignal.displ = s_filter(psignal.displ, dt, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)

    psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

    return psignal

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
# enf of get_bands

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
        print("[ERROR]: error loading her file. ")
        return False

    samples = dis_ns.size
    dt = time[1]

    # samples, dt, data, acceleration, velocity, displacement
    psignal_ns = seism_psignal(samples, dt, np.c_[dis_ns, vel_ns, acc_ns],
                               'c', acc_ns, vel_ns, dis_ns)
    psignal_ew = seism_psignal(samples, dt, np.c_[dis_ew, vel_ew, acc_ew],
                               'c', acc_ew, vel_ew, dis_ew)
    psignal_up = seism_psignal(samples, dt, np.c_[dis_up, vel_up, acc_up],
                               'c', acc_up, vel_up, dis_up)

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of read_file

def read_stamp(filename):
    """
    Get the time stamp from file's header
    """
    try:
        with open(filename) as input_file:
            try:
                header = input_file.readline().split()
                stamp = header[4].split(',')[-1].split(':')
                # tmp = stamp[2].split('.')
                # stamp[2] = tmp[0]
                # stamp.append(tmp[1])

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
# end of read_stamp

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
