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

    record_type = {'a': 'acceleration', 'v': 'velocity', 'd': 'displacement'}

    def set_dt(self, delta_t):
        if not isinstance(delta_t, float):
            print "\nError with time-step (dt).\n"
            return
        self.dt = delta_t

    def set_samples(self, samples):
        if not isinstance(samples, int):
            print "\nError with samples type.\n"
            return
        self.samples = samples

    def set_type(self, stype):
        if not isinstance(stype, str) or stype not in self.record_type:
            print "\nError with signal type (must be: a, v, or d).\n"
            return
        self.stype = stype

    def __init__(self, *args, **kwargs):
        """
        Initialize the record structure with all attributes empty.
        """
        self.stype = str()
        self.samples = int()
        self.dt = float()
        self.data = np.array([], float)

    def set_samples