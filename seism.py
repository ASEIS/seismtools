__author__ = 'rtaborda'

import numpy as np


class seism_record():
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

    def set_data(self, sdata):
        if not isinstance(sdata, np.array):
            print "\nError with signal data: not a numpy array.\n"
            return 3
        self.sdata = sdata
    #end set_data

    def set_type(self, stype):
        if not isinstance(stype, str) or stype not in self.record_type:
            print "\nError with signal type (must be: a, v, or d).\n"
        self.stype = stype
    #end set_type

    def __init__(self, *args, **kwargs):
        """
        Initialize the record structure with all attributes empty.
        If no parameters then all attributes set to None
        Can be initialize with unlabeled and labeld arguments
        Correct order for unlabeled arguments is: samples, dt, data, and signal type
        Acceptable
        """

        # Initialize to default values
        self.samples = None
        self.dt = None
        self.sdata = None
        self.stype = None

        if len(args) < 2:
            self.set_samples(args[0])
        if len(args) < 3:
            self.set_dt(args[1])
        if len(args) < 4:
            self.set_data(args[2])
        if len(args) < 5:
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
        print self.record_type[self.stype]
        print "Samples: "
        print "Delta t: %f" % self.dt
        print "Data:"
        print self.sdata
    #end __repr__
