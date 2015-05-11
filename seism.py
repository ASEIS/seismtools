#!/usr/bin/env python

__author__ = 'rtaborda'

import numpy as np
import matplotlib.pyplot as plt


class seism_signal():
    """
    This class implements attributes and methods related to
    a single seismic record. Here, a record is understood as
    a signal or seismogram, which can be any of three types:
    an acceleration, velocity, or displacement record with
    a given sampling rate size (delta t) and a total number
    of samples
    """

    record_type = {'a': 'Acceleration', 'v': 'Velocity', 'd': 'Displacement'}

    def set_samples(self, samples):
        if not isinstance(samples, int):
            print "\nError with samples type.\n"
        self.samples = samples
    #end set_samples

    def set_dt(self, delta_t):
        if not isinstance(delta_t, float):
            print "\nError with time-step (dt).\n"
        self.dt = delta_t
    #end set_dt

    def set_data(self, data):
        # check if the data passed is a numpy array
        # if not isinstance(sdata, np.array):
        if not isinstance(data, np.ndarray): 
            print "\nError with signal data: not a numpy array.\n"
        #     return 3
        self.data = data
    #end set_data

    def set_type(self, stype):
        if not isinstance(stype, str) or stype not in self.record_type:
            print "\nError with signal type (must be: a, v, or d).\n"
        self.type = stype
    #end set_type

    def __init__(self, *args, **kwargs):
        """
        Initialize the record structure with all attributes empty.
        Can be initialize with unlabeled and labeld arguments
        Correct order for unlabeled arguments is: samples, dt, data, and signal type
        If no parameters, samples set to 0, dt set to 0.01, type is accelerogram and array is empty
        Acceptable
        """

        # lines used for debugging
        # print 'args:   ', args
        # print 'kwargs: ', kwargs

        # Initialize to default values
        self.samples = 0
        self.dt = 0.01
        self.data = np.array([],float)
        self.type = 'a'

        if len(args) > 0:
            args_range = range(len(args))
            if 0 in args_range:
                self.set_samples(args[0])
            if 1 in args_range:
                self.set_dt(args[1])
            if 2 in args_range:
                self.set_data(args[2])
            if 3 in args_range:
                self.set_type(args[3])
                # all arguments were given in unlabled format
                return

        if len(kwargs) > 0:
            if 'samples' in kwargs:
                self.set_samples(kwargs['samples'])
            if 'dt' in kwargs:
                self.set_dt(kwargs['dt'])
            if 'data' in kwargs:
                self.set_data(kwargs['data'])
            if 'type' in kwargs:
                self.set_type(kwargs['type'])

        return
    #end __init__

    def __repr__(self):
        return "Seismic record: " + \
            self.record_type[self.type] + "\n" + \
            "Samples: %i" % self.samples + "\n" + \
            "Delta t: %.3f" % self.dt + "\n" + \
            "Data:\n" + '\n'.join(map(str, self.data))
    #end __repr__

    def plot(self, flag):
        t = np.arange(0, self.samples*self.dt, self.dt)
        plt.plot(t, self.data)
        if flag == 's':
            plt.show()
        elif flag == 'p':
            #TODO: print vector file
            pass
        return
    #end plot
#end signal class


class record(seism_signal):
    """
    This class extends the signal class to have addtitional
    attributes regarding time stamp and orientation
    """
    # TODO: to improve
    def __init__(self, station, location_lati, location_longi, hour, minute, seconds, fraction, tzone, orientation):
        self.set_tstamp(hour, minute, seconds, fraction, tzone)
        self.set_station(station, location_lati, location_longi)
        self.set_orientation(orientation)
    
    def set_station(self, station, location_lati, location_longi):
        # checking station name 
        if not isinstance(station, str): 
            print "\n**Error with station name.**\n"

        # checking latitude format being float+S/N
        if not isinstance(location_lati, str): 
            print "\n**Error with location latitude (Invalid instance type).**\n"
        elif location_lati[-1] not in ["N", "S"]: 
            print "\nError with location latitude (Invalid format).\n"
        else: 
            try:
                float(location_lati[0:-2])
            except ValueError:
                print "\nError with location latitude (Invalid format).\n"

        # checking longitude format being float+E/W
        if not isinstance(location_longi, str): 
            print "\n**Error with location longitude (Invalid instance type).**\n"
        elif location_longi[-1] not in ["E", "W"]: 
            print "\n**Error with location longitude (Invalid format).**\n"
        else: 
            try:
                float(location_longi[0:-2])
            except ValueError:
                print "\n**Error with location longitude (Invalid format).**\n"

        self.station = station 
        self.location_lati = location_lati 
        self.location_longi = location_longi
    #end set_station


    def set_orientation(self, orientation):
        # if the orientation is string, it should be either Up or Down 
        # if the orientation is int, it should between 0 and 360 
        if isinstance(orientation, str) and orientation in ["Up", "Down"]: 
            self.orientation = orientation
        elif isinstance(orientation, int) and orientation <= 360 and orientation >= 0:
            self.orientation = orientation
        else: 
            print "\n**Error with orientation (Invalid orientation).**\n"
    #end set_orientation

    def set_tstamp(self, hour, minute, seconds, fraction, tzone):
    	if not isinstance(hour, float):
    		print "\n**Error with record start time: hour.**\n"
        if not isinstance(minute, float):
            print "\n**Error with record start time: minute.**\n"
    	if not isinstance(seconds, float):
    		print "\n**Error with record start time: seconds.**\n"
    	if not isinstance(fraction, float):
    		print "\n**Error with record start time: fraction.**\n"
    	if not isinstance(tzone, str):
    		print "\n**Error with record start time: tzone.**\n"

        self.hour = hour
        self.minute = minute
        self.seconds = seconds 
        self.fraction = fraction 
        self.tzone = tzone
    # end set_tstamp

    # to test with record object
    def print_attr(self):
        print "station name: " + self.station
        print "station latitude: " + self.location_lati
        print "station longitude: " + self.location_longi
        print "tzone: " + self.tzone
        print "hour: " + str(self.hour) 
        print "minute: " + str(self.minute) 
        print "seconds: " + str(self.seconds) 
        print "fraction: " + str(self.fraction)
        print "depth: ??" 
        print "orientation: " + str(self.orientation)
# end record class
