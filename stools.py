#!/usr/bin/env python
# ===================================================================================
# The program contains general functions what may be used by other programs. 
# Including: filter; integral; derivative; FAS; 
# ===================================================================================
import numpy as np
import math 
from scipy.signal import filtfilt, ellip

fmin = 0.05 
fmax = 5.0 

def integrate(data, dt):
	data = np.cumsum(data*dt) #integrate
	data = ellip_filter(data, dt) #filt  
	# data = ellip_filter(np.cumsum(data*self.dt), self.dt, Wn = 0.075/((1.0/self.dt)/2.0), N = 7)
	return data


def derivative(data, dt):
	"""compute derivative of an numpy."""
	data = np.insert(data, 0, data[0])
	data = np.diff(data/dt)
	data = ellip_filter(data, dt) #filt
	return data 


def ellip_filter(data, dt, *args, **kwargs):
    """
    Correct order for unlabeled arguments is: data, dt, order N, rp, rs, Wn 
    """
    if not isinstance(data, np.ndarray): 
        print "\n[ERROR]: data is not an numpy array.\n"
        return 
    data = data 
    dt = dt
    N = 5
    rp = 0.1
    rs = 100 
    Wn = 0.05/((1.0/dt)/2.0)

    if len(args) > 0:
        args_range = range(len(args))
        if 0 in args_range:
            N = args[0]
        if 1 in args_range:
            rp = args[1]
        if 2 in args_range:
            rs = args[2]
        if 3 in args_range:
            Wn = args[3]

    if len(kwargs) > 0:
        if 'N' in kwargs:
            N = kwargs['N']
        if 'rp' in kwargs:
            rp = kwargs['rp']
        if 'rs' in kwargs:
            rs = kwargs['rs']
        if 'Wn' in kwargs:
            Wn = kwargs['Wn']

    # create highpass elliptic filter 
    b, a = ellip(N = N, rp = rp, rs = rs, Wn = Wn, btype = 'highpass', analog=False)
    data = filtfilt(b, a, data)
    return data 

def smooth(data, factor): 
    # factor = 3; c = 0.5, 0.25, 0.25
    # TODO: fix coefficients for factors other than 3 
    c = 0.5/(factor-1)
    for i in range(1, data.size-1):
        data[i] = 0.5*data[i] + c*data[i-1] + c*data[i+1]
    return data 


def FAS(data, dt, points, fmin, fmax, s_factor):
    afs = abs(np.fft.fft(data, points))*dt
	# freq = (1/signal.dt)*range(points)/points
    freq = (1/dt)*np.array(range(points))/points

    deltaf = (1/dt)/points

    inif = int(fmin/deltaf)
    endf = int(fmax/deltaf) + 1
    
    afs = afs[inif:endf]
    afs = smooth(afs, s_factor)
    freq = freq[inif:endf]
    return freq, afs

    # afs = afs[inif:endf]
    # afs = smooth(afs, s_factor) # smooth the curve 
    # freq = freq[inif:endf]
    # return freq, afs 

def get_points(samples1, samples2):
	# points is the least base-2 number that is greater than max samples 
	power = int(math.log(max(samples1, samples2), 2)) + 1 
	return 2**power 
# end of get_points

def set_bound(f1, f2): 
	"""to get fmin and fmax from other programs"""
	global fmin 
	global fmax 
	fmin = f1 
	fmax = f2

	print fmin, fmax 
# end of set_bound

