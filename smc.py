#!/usr/bin/env python
from __future__ import division
from seism import *
from scipy.signal import filtfilt, ellip
import os
import math 

discard = {'dam': 'Dam', 'Fire Sta': 'Fire Station', 'Acosta Res': 'Acosta Res', 'Bldg': 'Building', 'Br': 'Interchange Bridge'}
destination = '' 
def get_destination(d):
    """The function is to get the user input from process.py."""
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
    for i in range(1,len(channels)):
        del channels[i][0]

    for i in range(len(channels)):

        # check this is the uncorrected acceleration data
        ctype = channels[i][0][0:24].lower()
        if ctype != "uncorrected accelerogram":
            return channels, 0
        else:
            dtype = 'a'


        # tmp = channels[i][3].split('.')
        # network = tmp[1]
        # station_id = tmp[2]
        network = filename.split('/')[-1].split('.')[0][0:2].upper()
        station_id = filename.split('/')[-1].split('.')[0][2:].upper()


        # get location's latitude and longitude 
        tmp = channels[i][4].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        depth = 0.0


        # get station name
        station = channels[i][5][0:40].strip()

        # get data type 
        # tmp = channels[i][5].split()
        # dtype = tmp[-1]
        # if dtype == "Acceleration": 
        # 	dtype = 'a'
        # elif dtype == "Velocity": 
        # 	dtype = 'v'
        # elif dtype == "Displacement": 
        # 	dtype = 'd'
        # else: 
        # 	dtype = "Unknown"


        # get orientation, convert to int if it's digit 
        tmp = channels[i][6].split()
        orientation = tmp[2]
        if orientation.isdigit():
            orientation = int(orientation)
        location = channels[i][6][36:].strip()
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
        time = hour + ":" + minute + ":" + seconds + "." + fraction + " " + tzone
   

        # get number of samples and dt 
        tmp = channels[i][27].split()
        samples = int(tmp[0])
        dt = 1/int(tmp[4])

        # get signals' data 
        tmp = channels[i][28:]
        signal = str()
        for s in tmp:
            signal += s
        data = process_signal(signal)


        # record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth = depth, 
            orientation = orientation, date = date, time = time)
        # record.print_attr()

        record_list.append(record)

    # return a list of records, and corresponding network code and station id 
    return record_list, network, station_id

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
    # this row corresponds to the channel delimiter
    for i in range(1,len(channels)):
        del channels[i][0]

    for i in range(len(channels)):

        tmp = channels[i][0].split()

        # check this is the uncorrected acceleration data
        ctype = (tmp[0] + " " + tmp[1]).lower()
        if ctype != "corrected accelerogram":
            return channels, 0

        # get network code and station id 
        # network = tmp[2].split('.')[1]
        # station_id = tmp[2].split('.')[2]
        network = filename.split('/')[-1].split('.')[0][0:2].upper()
        station_id = filename.split('/')[-1].split('.')[0][2:].upper()

        # get orientation, convert to int if it's digit 
        # orientation = tmp[5]
        # if orientation.isdigit():
        #     orientation = int(orientation)


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
        location = channels[i][7][36:].strip()
        print location 
        # if 'Depth' in location: 
        #     depth = float(location.split()[2])
        # else:
        #     pass 
            # TODO: set location 


        # get station name
        station = channels[i][6][0:40].strip()


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
        time = hour + ":" + minute + ":" + seconds + "." + fraction + " " + tzone
   

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

        a_data = process_signal(a_signal)
        v_data = process_signal(v_signal)
        d_data = process_signal(d_signal)
        data = np.c_[a_data, v_data, d_data]

        record = seism_record(samples, dt, a_data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        precord = seism_precord(samples, dt, data, dtype, accel = a_data, displ = d_data, velo = v_data, orientation = orientation, date = date, time = time, depth = depth,
            latitude = location_lati, longitude = location_longi)
        if orientation in [0, 360, 180, -180]:
            orientation = 'N'
        elif orientation in [90, 270, -90, -270]:
            orientation = 'E'
        elif orientation in ['Up', 'Down']:
            orientation = 'Z'
        else:
            orientation = ' '

        print_smc(network + "." + station_id + ".V2" + orientation + ".txt", record)
        record_list.append(precord)

    return record_list, network, station_id



def process_signal(signal):
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


def print_smc(filename, record):
    """
    The function generates .txt files for each channel/record 
    """
    global destination
    network = filename.split('.')[0]
    station = filename.split('.')[1]
    info = "XXX"  # info = band + instrument type + orientation 
    
    # generate a text file (header + data)
    header = "# " + network + " " + station + " " + info + " " + record.date + "," + record.time + " " + str(record.samples) + " " + str(record.dt) + "\n"
    try:
        f = open(destination + '/' + filename, 'w')
    except IOError, e:
        print e
        # return 
    
    f.write(header)
    # descriptor = '{:>12.7f}' + '\n'
    descriptor = '{:>f}' + '\n'
    if record.data.size != 0: 
        for d in np.nditer(record.data):
            # f.write(str(d)+"\n")
            f.write(descriptor.format(float(d)))
    f.close()
    print "*Generated .txt file at: " + destination + "/" + filename

#end of print_smc

def print_her(filename, record_list):
    """
    The function generates .her files for each station (with all three channels included)
    """
    global destination
     # if there are more than three channels, save for later 
    if len(record_list) > 3:
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

    # round data to 7 decimals in order to print properly 
    for precord in record_list:
        if precord.orientation == " ":
            print "[ERROR]: missing orientation"
            orientation = " "
        elif precord.orientation in [0, 360, 180, -180]:
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
             # handling degrees such as 60, 120 etc. 
             pass

    # get a list of time incremented by dt 
    time = [0.000]
    samples = precord.samples 
    while samples > 1:
        time.append(time[len(time)-1] + precord.dt)
        samples -= 1 
        descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    f.write(descriptor.format("time", "dis_ns", "dis_ew", "dis_up", "vel_ns", "vel_ew", "vel_up", "acc_ns", "acc_ew", "acc_up")) # header 
    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew, vel_up, acc_ns, acc_ew, acc_up):
        f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 ))
    f.close()
    print "*Generated .her file at: " + destination + "/" + filename
#end of print_her 

def rotate(record_list):
    """
    The function is to transfrom data for channels with special orientations. 
    """
    tmp = []
    for record in record_list:
        if record.orientation != 'Up':
            tmp.append(record)
            record_list.remove(record)
    if len(tmp) != 2: 
        return False 
    x = tmp[0].orientation
    y = tmp[1].orientation
    if x > y: 
        list(reversed(tmp))

    # rotate 
    matrix = np.array([(math.cos(math.radians(x)), -math.sin(math.radians(x))), (math.sin(math.radians(x)), math.cos(math.radians(x)))])
    data = matrix.dot([tmp[0].data, tmp[1].data])

    # transform the first record with North orientation 
    tmp[0].data = data[0]
    tmp[0].orientation = 0

    # transform the second record with East orientation
    tmp[1].data = data[1]
    tmp[1].orientation = 90

    record_list += tmp
    return record_list
# end of rotate

def print_warning(warning):
    """
    The function is to generate warning messages for channles at special locations. 
    """
    f = open(destination + '/warning.txt', 'a')
    f.write(warning +"\n")
    f.close()
# end of print_warning

def process_record_list(network, station_id, record_list):
    """
    The function is to take a list of V1 [record]s, then use their data to get velocity and displacement,
    then create [precord] objects, finally return a list of them. 
    If encounter file with more than three channels, print message and return False.
    If encounter special orientations, rotate by matrix, and recursivly calling the function to continue processing. 
    """
    global discard
    rotate_flag = False 
    if len(record_list) != 3:
        print "==[The function is processing files with 3 channels only.]=="
        return False 
    processed_list = []

    for record in record_list:
        # process data of record 
        filename = record.process_smc_v1(network, station_id)
        if filename == False: # if encounter special orientations. 
            rotate_flag = True 
            break 

        # generate .txt file for record 
        print_smc(filename, record)

        # get velocity and displacement
        velocity = record.integrate(record.data)
        displacement = record.integrate(velocity)
        precord = seism_precord(record.samples, record.dt, record.data, record.type, accel = record.data, displ = displacement, velo = velocity, 
            orientation = record.orientation, date = record.date, time = record.time, depth = record.depth, latitude = record.location_lati, longitude = record.location_longi)
        processed_list.append(precord)

    # rotation 
    if rotate_flag:
        record_list = rotate(record_list)
        if record_list == False: 
            return False 
        return process_record_list(network, station_id, record_list) # recursively calling the function to continue processing

    # if the channel is in the location should be discarded; print warning 
    for key in discard:
            if key in record_list[0].station:
                filename = network + station_id + '.V1'
                warning = filename + " was processed, but it's from " + discard[key]
                print_warning(warning)
                break 

    return processed_list
#end process_record_list 


# =================================== Test Cases ====================================
# record_list, network, station_id = load_smc_v1('CIRSB.D2.RAW')
# filename = network + "." + station_id + ".V1.her"
# print_her(filename, process_record_list(network, station_id, record_list))


# record_list, network, station_id = load_smc_v2('CE14783.V2')
# filename = network + "." + station_id + ".V2.her"
# print_her(filename, record_list)


# record_list, network, station_id = load_smc_v1('NP05396.RAW')
# filename = network + "." + station_id + ".V1.her"
# print_her(filename, process_record_list(network, station_id, record_list))







