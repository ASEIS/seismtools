#!/usr/bin/env python
from __future__ import division
from seism import *
from scipy.signal import filtfilt, ellip

def load_smc_v1(filename):
    record_list = []
    
    # loads station into a string
    fp = open(filename, 'r')
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


        tmp = channels[i][3].split('.')
        network = tmp[1]
        station_id = tmp[2]

        # get location's latitude and longitude 
        tmp = channels[i][4].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        depth = 0.0


        # get station name
        station = channels[i][5][0:40].strip()

        # get data type 
        tmp = channels[i][5].split()
        dtype = tmp[-1]
        if dtype == "Acceleration": 
        	dtype = 'a'
        elif dtype == "Velocity": 
        	dtype = 'v'
        elif dtype == "Displacement": 
        	dtype = 'd'
        else: 
        	dtype = "Unknown"


        # get orientation, convert to int if it's digit 
        tmp = channels[i][6].split()
        orientation = tmp[2]
        if orientation.isdigit():
            orientation = int(orientation)

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

    # return a list of records and corresponding network code and station id 
    return record_list, network, station_id


def process_smc_v1(record_list, network, station_id):
    """
    The function process a list of records by converting orientation and multiplying data 
    """
    # if there are more than three channels, save for later 
    if len(record_list) > 3:
        print "==[The function is processing files with 3 channels only.]=="
        return 

    else: 
        for record in record_list:
            # ======================================= processing orientation ========================================== 
            # If the orientation was not set properly, it would be empty string by default 
            if record.orientation == " ":
                print "[ERROR]: missing orientation"
                orientation = " "
            elif record.orientation in [0, 360]:
                orientation = 'N'
            elif record.orientation in [180, -180]:
                orientation = 'N'
                record.data = record.data*-1
            elif record.orientation in [90, -270]:
                orientation = 'E'
            elif record.orientation in [-90, 270]:
                orientation = 'E'
                record.orientation = record.data*-1
            elif record.orientation == "Up" or record.orientation == "Down":
                orientation = 'Z'
            else: 
                # handling degrees such as 60, 120 etc. 
                pass

            # ======================================= processing data ================================================  
            # create highpass elliptic filter 
            b, a = ellip(N = 5, rp = 0.1, rs = 100, Wn = 0.05/((1.0/record.dt)/2.0), btype = 'highpass', analog=False)
            filtered_signal = filtfilt(b, a, record.data)
            filtered_record = seism_record(record.samples, record.dt, filtered_signal, record.type, record.station, record.location_lati, record.location_longi, depth = record.depth, 
            orientation = record.orientation, date = record.date, time = record.time)

            # record.plot('s')
            # filtered_record.plot('s')

            # minus average and multiply by 981 
            filtered_record.data = 981*(filtered_record.data - np.average(filtered_record.data))

            filename = network + "." + station_id + ".V1" + orientation + ".txt"
            print_smc(filename, filtered_record)
        return 
    return


def print_smc(filename, record):
    """
    The function generate files for each channel/record 
    """
    # generate a text file (header + data)
    header = "# " + record.date + " " + record.time + " Samples: " + str(record.samples) + " dt: " + str(record.dt) + "\n"
    f = open(filename, 'w')
    f.write(header)
    for d in np.nditer(record.data):
        f.write(str(d)+"\n")
    f.close()
#end of print_smc_v1()



def load_smc_v2(filename):
    record_list = []
    
    # loads station into a string
    fp = open(filename, 'r')
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
        network = tmp[2].split('.')[1]
        station_id = tmp[2].split('.')[2]

        # get orientation, convert to int if it's digit 
        orientation = tmp[5]
        if orientation.isdigit():
            orientation = int(orientation)


        # get location's latitude and longitude 
        tmp = channels[i][5].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        depth = 0.0


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

        v_data = process_signal(a_signal)
        a_data = process_signal(v_signal)
        d_data = process_signal(d_signal)

        # Objects for current testing 
        record1 = seism_record(samples, dt, a_data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        # record2 = seism_record(samples, dt, v_data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        # record3 = seism_record(samples, dt, d_data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        if orientation in [0, 360, 180, -180]:
            orientation = 'N'
        elif orientation in [90, 270, -90, -270]:
            orientation = 'E'
        elif orientation in ['Up', 'Down']:
            orientation = 'Z'
        else:
            orientation = ' '
        filename = network + "." + station_id + ".V2" + orientation + ".txt"
        print_smc(filename, record1)


def process_signal(signal):
    """
    The function is to convert signal into an numpy array of float numbers 
    """
    # avoid negative number being stuck  
    signal = signal.replace('-', ' -')
    signal = signal.split()

    data = []
    for s in signal:
        data.append(float(s))
    data = np.array(data)
    return data 


load_smc_v2('NCNHC.V2')

# test with two files 
record_list, network, station_id = load_smc_v1('NCNHC.V1')
process_smc_v1(record_list, network, station_id)

# record_list, network, station_id = load_smc_v1('CIQ0028.V1')
# process_smc_v1(record_list, network, station_id)


# TODO:  
# 3. three column array 
# seism-psignal: sub of signal 
# seism-precord: check station, orientation etc; sub of record 
