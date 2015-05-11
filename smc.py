#!/usr/bin/env python
from __future__ import division
from seism import *

def load_smc_v1(filename):
    
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


        # get location's latitude and longitude 
        tmp = channels[i][4].split()
        location_lati = tmp[3][:-1]
        location_longi = tmp[4]
        # depth = ? 


        # get station name
        station = channels[i][5][0:40]

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


        # get orientation 
        tmp = channels[i][6].split()
        orientation = tmp[2]

        # get time 
        # start_time = channels[i][3][37:80]
        tmp = channels[i][14].split()
        hour = float(tmp[0])
        minute = float(tmp[1])
        seconds = float(tmp[2])
        fraction = float(tmp[4])
        tzone = channels[i][3].split()[-2]
   

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



        # record = seism_signal(samples=samples, dt=dt, data=data, type=dtype)
        # record.plot('s')
        # print record

        print "====================================================="
        print "channel: " + ctype
        print "samples: " + str(samples)
        print "dt: " + str(dt) 
        print "data type: " + dtype 
        print "station name: " + station
        print "station latitude: " + location_lati
        print "station longitude: " + location_longi
        print "tzone: " + tzone
        print "hour: " + str(hour) 
        print "minute: " + str(minute) 
        print "seconds: " + str(seconds) 
        print "fraction: " + str(fraction)
        print "depth: ??" 
        print "orientation: " + orientation

        # print data
        pass

    return channels, 1

# test with two files 
load_smc_v1('CIQ0028.V1')
load_smc_v1('NCNHC.V1')

# TODO: 
# 1. check fraction second 
# 4. depth
# 5. convert orientation degree to directions? or convert up/down to int? 
# 6. data? 

