#!/usr/bin/env python
from __future__ import division
import os
import numpy as np
from seism import seism_record, seism_station, seism_precord

discard = {'dam': 'Dam', 'Fire Sta': 'Fire Station',
           'Acosta Res': 'Acosta Res', 'Bldg': 'Building',
           'Br': 'Interchange Bridge'}
destination = ''

def set_destination(d):
    """
    The function is to set the destination variable with the user input.
    """
    global destination
    destination = d

def load_smc_v1(filename):
    record_list = []

    # loads station into a string
    try:
        fp = open(filename, 'r')
    except IOError, e:
        print e
        # return

    channels = fp.read()
    fp.close()

    # splits the string by channels
    channels = channels.split('/&')
    del(channels[len(channels)-1])

    # splits the channels
    for i in range(len(channels)):
        channels[i] = channels[i].split('\r\n')

    # clean the first row in all but the first channel
    # this row corresponds to the channel delimiter
    for i in range(1, len(channels)):
        del channels[i][0]

    for i in range(len(channels)):

        # check this is the uncorrected acceleration data
        ctype = channels[i][0][0:24].lower()
        if ctype != "uncorrected accelerogram":
            print "[ERROR]: processing uncorrected accelerogram ONLY."
            return False
        else:
            dtype = 'a'

        network = filename.split('/')[-1].split('.')[0][0:2].upper()
        station_id = filename.split('/')[-1].split('.')[0][2:].upper()

        # get location's latitude and longitude
        tmp = channels[i][4].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        depth = 0.0

        # get station name
        station_name = channels[i][5][0:40].strip()

        # get orientation, convert to int if it's digit
        tmp = channels[i][6].split()
        orientation = tmp[2]
        if orientation.isdigit():
            orientation = int(orientation)
        # location = channels[i][6][36:].strip()
        # if 'Depth' in location:
        #     depth = float(location.split()[2])
        # else:
        #     pass
            # TODO: set location

        # get date and time; set to fixed format
        start_time = channels[i][3][37:80].split()
        date = start_time[2][:-1]

        tmp = channels[i][14].split()
        hour = tmp[0]
        minute = tmp[1]
        seconds = tmp[2]
        # fraction = tmp[4]
        fraction = tmp[3]
        tzone = channels[i][3].split()[-2]
        time = "%s:%s:%s.%s %s" % (hour, minute, seconds, fraction, tzone)

        # get number of samples and dt
        tmp = channels[i][27].split()
        samples = int(tmp[0])
        dt = 1/int(tmp[4])

        # get signals' data
        tmp = channels[i][28:]
        signal = str()
        for s in tmp:
            signal += s
        data = read_data(signal)

        record = seism_record(samples, dt, data, dtype, station_name,
                              location_lati, location_longi, depth=depth,
                              orientation=orientation,
                              date=date, time=time)

        record_list.append(record)

    station = seism_station(record_list, network, station_id, 'V1')
    return station

def load_smc_v2(filename):
    record_list = []

    # loads station into a string
    try:
        fp = open(filename, 'r')
    except IOError, e:
        print e
        # return

    channels = fp.read()
    fp.close()

    # splits the string by channels
    channels = channels.split('/&')
    del(channels[len(channels)-1])

    # splits the channels
    for i in range(len(channels)):
        channels[i] = channels[i].split('\r\n')

    # clean the first row in all but the first channel
    for i in range(1, len(channels)):
        del channels[i][0]

    for i in range(len(channels)):

        tmp = channels[i][0].split()

        # check this is the uncorrected acceleration data
        ctype = (tmp[0] + " " + tmp[1]).lower()
        if ctype != "corrected accelerogram":
            print "[ERROR]: processing corrected accelerogram ONLY."
            return False

        # get network code and station id
        network = filename.split('/')[-1].split('.')[0][0:2].upper()
        station_id = filename.split('/')[-1].split('.')[0][2:].upper()

        # get location's latitude and longitude
        tmp = channels[i][5].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        depth = 0.0

        # get orientation, convert to int if it's digit
        tmp = channels[i][7].split()
        orientation = tmp[2]
        if orientation.isdigit():
            orientation = int(orientation)
        # location = channels[i][7][36:].strip()

        # if 'Depth' in location:
        #     depth = float(location.split()[2])
        # else:
        #     pass
            # TODO: set location

        # get station name
        station_name = channels[i][6][0:40].strip()

        # get date and time; set to fixed format
        start_time = channels[i][4][37:80].split()
        date = start_time[2][:-1]

        tmp = channels[i][26].split()
        hour = tmp[0]
        minute = tmp[1]
        seconds = tmp[2]
        # fraction = tmp[4]
        fraction = tmp[3]
        tzone = channels[i][4].split()[-2]
        time = "%s:%s:%s.%s %s" % (hour, minute, seconds, fraction, tzone)

        # get number of samples and dt
        tmp = channels[i][45].split()
        samples = int(tmp[0])
        dt = float(tmp[8])

        # get signals' data
        tmp = channels[i][45:]
        a_signal = str()
        v_signal = str()
        d_signal = str()
        for s in tmp:
            # detecting separate line and get data type
            if "points" in s:
                line = s.split()
                if line[3] == "accel":
                    dtype = 'a'
                elif line[3] == "veloc":
                    dtype = 'v'
                elif line[3] == "displ":
                    dtype = 'd'
                else:
                    dtype = "Unknown"

            # processing data
            else:
                if dtype == 'a':
                    a_signal += s
                elif dtype == 'v':
                    v_signal += s
                elif dtype == 'd':
                    d_signal += s

        a_data = read_data(a_signal)
        v_data = read_data(v_signal)
        d_data = read_data(d_signal)
        data = np.c_[d_data, v_data, a_data]

        precord = seism_precord(samples, dt, data, 'c', station_name,
                                accel=a_data, displ=d_data, velo=v_data,
                                orientation=orientation, date=date,
                                time=time, depth=depth,
                                latitude=location_lati,
                                longitude=location_longi)
        record_list.append(precord)

    station = seism_station(record_list, network, station_id, 'V2')
    if not station.list:
        return False
    else:
        return station
# end of load_smc_v2

def read_data(signal):
    """
    The function is to convert signal data into an numpy array of float numbers
    """
    # avoid negative number being stuck
    signal = signal.replace('-', ' -')
    signal = signal.split()

    data = []
    for s in signal:
        data.append(float(s))
    data = np.array(data)
    return data

def print_smc(station):
    """
    The function generates .txt files for each channel/record
    """

    orientation = ''

    for record in station.list:
        if record.orientation in [0, 180, 360, -180]:
            orientation = 'N'
        elif record.orientation in [90, -90, 270, -270]:
            orientation = 'E'
        elif record.orientation in ['Up', 'Down']:
            orientation = 'Z'

        filename = "%s.%s.%s.txt" % (station.network, station.id,
                                     station.type + orientation)

        # generate a text file (header + data)
        header = "# %s %s %s %s,%s %s %s\n" % (station.network, station.id,
                                               station.type + orientation,
                                               record.date, record.time,
                                               str(record.samples),
                                               str(record.dt))
        try:
            f = open(os.path.join(destination, filename), 'w')
        except IOError, e:
            print e
            # return

        f.write(header)
        descriptor = '{:>f}' + '\n'

        if record.accel.size != 0:
            for d in np.nditer(record.accel):
                f.write(descriptor.format(float(d)))
        f.close()
        print "*Generated .txt file at: " + os.path.join(destination, filename)
#end of print_smc

def print_bbp(station):
    """
    This function generates .bbp files for
    each of velocity/acceleration/displacement
    """

    filename_base = "%s.%s.%s" % (station.network, station.id, station.type)

    # round data to 7 decimals in order to print properly
    for precord in station.list:
        if precord.orientation in [0, 360, 180, -180]:
            dis_ns = precord.displ.tolist()
            vel_ns = precord.velo.tolist()
            acc_ns = precord.accel.tolist()
        elif precord.orientation in [90, -270, -90, 270]:
            dis_ew = precord.displ.tolist()
            vel_ew = precord.velo.tolist()
            acc_ew = precord.accel.tolist()
        elif precord.orientation == "Up" or precord.orientation == "Down":
            dis_up = precord.displ.tolist()
            vel_up = precord.velo.tolist()
            acc_up = precord.accel.tolist()
        else:
            pass

    # Prepare to output
    out_data = [['dis', dis_ns, dis_ew, dis_up, 'displacement', 'm'],
                ['vel', vel_ns, vel_ew, vel_up, 'velocity', 'm/s'],
                ['acc', acc_ns, acc_ew, acc_up, 'acceleration', 'm/s^2']]

    for data in out_data:
        filename = "%s.%s.bbp" % (filename_base, data[0])
        try:
            out_fp = open(os.path.join(destination, filename), 'w')
        except IOError, e:
            print e
            continue

        # Start with time = 0.0
        time = [0.000]
        samples = precord.samples
        while samples > 1:
            time.append(time[len(time)-1] + precord.dt)
            samples -= 1

        # Write header
        out_fp.write("# Station: %s_%s\n" % (station.network, station.id))
        out_fp.write("#    time= %s,%s\n" % (precord.date, precord.time))
        out_fp.write("#     lon= %s\n" % (station.longitude))
        out_fp.write("#     lat= %s\n" % (station.latitude))
        out_fp.write("#   units= %s\n" % (data[5]))
        out_fp.write("#\n")
        out_fp.write("# Data fields are TAB-separated\n")
        out_fp.write("# Column 1: Time (s)\n")
        out_fp.write("# Column 2: N/S component ground "
                     "%s (%s) (+ is 000)\n" % (data[4], data[5]))
        out_fp.write("# Column 3: E/W component ground "
                     "%s (%s) (+ is 090)\n" % (data[4], data[5]))
        out_fp.write("# Column 4: U/D component ground "
                     "%s (%s) (+ is upward)\n" % (data[4], data[5]))
        out_fp.write("#\n")

        # Write timeseries
        for val_time, val_ns, val_ew, val_ud in zip(time, data[1],
                                                    data[2], data[3]):
            out_fp.write("%5.7f   %5.9e   %5.9e    %5.9e\n" %
                         (val_time, val_ns, val_ew, val_ud))

        # All done, close file
        out_fp.close()
        print "*Generated .bbp file at: " + os.path.join(destination, filename)

def print_her(station):
    """
    The function generates .her files for each
    station (with all three channels included)
    """

    filename = '%s.%s.%s.her' % (station.network, station.id, station.type)

    try:
        f = open(os.path.join(destination, filename), 'w')
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
    # orientation = ''

    # round data to 7 decimals in order to print properly
    for precord in station.list:
        if precord.orientation in [0, 360, 180, -180]:
            dis_ns = precord.displ.tolist()
            vel_ns = precord.velo.tolist()
            acc_ns = precord.accel.tolist()
        elif precord.orientation in [90, -270, -90, 270]:
            dis_ew = precord.displ.tolist()
            vel_ew = precord.velo.tolist()
            acc_ew = precord.accel.tolist()
        elif precord.orientation == "Up" or precord.orientation == "Down":
            dis_up = precord.displ.tolist()
            vel_up = precord.velo.tolist()
            acc_up = precord.accel.tolist()
        else:
            pass
        # orientation = precord.orientation

    # get a list of time incremented by dt
    time = [0.000]
    samples = precord.samples
    while samples > 1:
        time.append(time[len(time)-1] + precord.dt)
        samples -= 1

    header = "# %s %s %s %s,%s %s %s" % (station.network, station.id,
                                         station.type, precord.date,
                                         precord.time, str(precord.samples),
                                         str(precord.dt))
    f.write(header)

    descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    f.write(descriptor.format("# time",
                              "dis_ns", "dis_ew", "dis_up",
                              "vel_ns", "vel_ew", "vel_up",
                              "acc_ns", "acc_ew", "acc_up")) # header

    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time,
                                                      dis_ns, dis_ew, dis_up,
                                                      vel_ns, vel_ew, vel_up,
                                                      acc_ns, acc_ew, acc_up):
        f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9))
    f.close()
    print "*Generated .her file at: " + os.path.join(destination, filename)
#end of print_her
