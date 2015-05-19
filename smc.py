#!/usr/bin/env python
from __future__ import division
from seism import *
from scipy.io import loadmat
from scipy.signal import butter, filtfilt, ellip
from matplotlib.pyplot import plot
from scipy import signal

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
        signal = signal.split()
         
         # make the signal a numpy array of float numbers
        data = []
        for s in signal: 
        	data.append(float(s)) 
        # print data 
        data = np.array(data)


        # record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth = depth, 
            orientation = orientation, date = date, time = time)
        # record.print_attr()
        # signal = seism_signal(dt=dt, samples=samples, data=data, type=dtype)
        # signal = seism_signal(samples,dt,data,dtype)
        # signal = seism_signal(samples,dt,type=dtype,data=data)

        # signal.plot('s')
        # record.plot('s')
        record_list.append(record)

    # return channels, 1
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
            # filtering data 
            # filter type = elliptic 
            # high pass filter 
            # 0.05Hz 
            # zero phrase 
            # matlab function filtfilt 


            print record.data.size
            b, a = ellip(N = 17, rp = 0.01, rs = 60, Wn = 0.05/((1.0/record.dt)/2.0), btype = 'highpass', analog=False)
            # b, a = ellip(N = 1, rp = 4, rs = 6, Wn = 0.075, btype = 'highpass', analog=False)
            filtered_signal = filtfilt(b, a, record.data)
            print filtered_signal.size 

            # minus average and multiply by 981 
            record.data = 981*(record.data - np.average(record.data))
            # record.plot('s')

            filename = network + "." + station_id + "." + orientation + ".txt"
            print_smc_v1(filename, record)
        return 
    return


def print_smc_v1(filename, record):
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




# test with two files 
record_list, network, station_id = load_smc_v1('NCNHC.V1')
process_smc_v1(record_list, network, station_id)

record_list, network, station_id = load_smc_v1('CIQ0028.V1')
process_smc_v1(record_list, network, station_id)


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

        # get data type 
        tmp = channels[i][6].split()
        dtype = tmp[-1]
        if dtype == "Acceleration": 
            dtype = 'a'
        elif dtype == "Velocity": 
            dtype = 'v'
        elif dtype == "Displacement": 
            dtype = 'd'
        else: 
            dtype = "Unknown"


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
        tmp = channels[i][46:]
        signal = str()
        for s in tmp:
            # excluding separate line 
            if not "points" in s: 
                signal += s

        # avoid negative number being stacked 
        signal = signal.replace('-', ' -')
        signal = signal.split()
         
         # make the signal a numpy array of float numbers
        data = []
        for s in signal: 
            data.append(float(s))
        data = np.array(data)


        # for d in np.nditer(data):
        #     print d 
            # pass
        # print data 


        # record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth, date, time, orientation)
        record = seism_record(samples, dt, data, dtype, station, location_lati, location_longi, depth = depth, 
            orientation = orientation, date = date, time = time)
        record.print_attr()

        # record.plot('s')
        record_list.append(record)

    # return a list of records and corresponding network code and station id 
    return record_list, network, station_id

# load_smc_v2('NCNHC.V2')

# TODO: 
# 1. dt 
# 2. error with plot 
# 3. 