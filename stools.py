#!/usr/bin/env python
# ===================================================================================
# The program contains general functions what may be used by other programs. 
# Including: filter; integral; derivative; FAS; 
# ===================================================================================
from __future__ import division
import numpy as np
import math 
from scipy.signal import filtfilt, ellip, butter, kaiser
from scipy.fftpack import fft, fftshift

def integrate(data, dt):
	data = np.cumsum(data*dt) 
	return data


def derivative(data, dt):
	"""compute derivative of an numpy."""
	data = np.insert(data, 0, data[0])
	data = np.diff(data/dt)
	return data 

def s_filter(*args, **kwargs):
    """
    correct order for unlabeled arguments is data, dt; 
    """
    data = np.array([],float)
    dt = 0.0 
    fami = {'ellip': ellip, 'butter': butter}

    if len(args) == 2:
        data = args[0]
        dt = args[1]
    else: 
        print "[ERROR]: filter missing data and dt."
        return data 

    if not isinstance(data, np.ndarray):
        print "[ERROR]: data input for filter is not an numpy array."
        return data 

    # default values 
    N = 5
    rp = 0.1
    rs = 100 
    Wn = 0.05/((1.0/dt)/2.0)
    fmin = 0.0 
    fmax = 0.0
    a = np.array([],float)
    b = np.array([],float)

    if len(kwargs) > 0: 
        if 'type' in kwargs:
            btype = kwargs['type']
        if 'N' in kwargs:
            N = kwargs['N']
        if 'rp' in kwargs:
            rp = kwargs['rp']
        if 'rs' in kwargs:
            rs = kwargs['rs']
        if 'Wn' in kwargs:
            Wn = kwargs['Wn']
        if 'fmin' in kwargs :
            fmin = kwargs['fmin']
            w_min = fmin/((1.0/dt)/2.0)
        if 'fmax' in kwargs:
            fmax = kwargs['fmax']
            w_max = fmax/((1.0/dt)/2.0)

        if fmin and fmax: 
            Wn = [w_min, w_max]
        elif fmax:
            Wn = w_max

        # calling filter 
        if kwargs['family'] == 'ellip': 
            b, a = fami[kwargs['family']](N = N, rp = rp, rs = rs, Wn = Wn, btype = btype, analog=False)
        elif kwargs['family'] == 'butter':
            b, a = fami[kwargs['family']](N = N, Wn = Wn, btype = btype, analog=False)

        data = filtfilt(b, a, data)
        return data 
# end of s_filter 


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


def get_points(samples1, samples2):
	# points is the least base-2 number that is greater than max samples 
	power = int(math.log(max(samples1, samples2), 2)) + 1 
	return 2**power 
# end of get_points

def get_period(tmin, tmax):
    """ Return an array of period T """
    # tmin = 1/fmax 
    # tmax = 1/fmin 
    a = np.log10(tmin)
    b = np.log10(tmax) 

    period = np.linspace(a, b, 20)
    period = np.power(10, period)
    return period 


def max_osc_response(acc, dt, csi, period, ini_disp, ini_vel):
    signal_size = acc.size 

    # initialize numpy arrays
    d = np.empty((signal_size))
    v = np.empty((signal_size))
    aa = np.empty((signal_size)) 

    d[0] = ini_disp
    v[0] = ini_vel

    w = 2*math.pi/period
    ww = w**2 
    csicsi = csi**2 
    dcsiw=2*csi*w

    rcsi=math.sqrt(1-csicsi)
    csircs=csi/rcsi
    wd=w*rcsi
    ueskdt=-1/(ww*dt)
    dcsiew=2*csi/w
    um2csi=(1-2*csicsi)/wd
    e=math.exp(-w*dt*csi)
    s=math.sin(wd*dt)
    c0=math.cos(wd*dt);
    aa[0]=-ww*d[0]-dcsiw*v[0]

    ca=e*(csircs*s+c0)
    cb=e*s/wd
    cc=(e*((um2csi-csircs*dt)*s-(dcsiew+dt)*c0)+dcsiew)*ueskdt
    cd=(e*(-um2csi*s+dcsiew*c0)+dt-dcsiew)*ueskdt
    cap=-cb*ww
    cbp=e*(c0-csircs*s)
    ccp=(e*((w*dt/rcsi+csircs)*s+c0)-1)*ueskdt
    cdp=(1-ca)*ueskdt

    for i in range(1, signal_size):
        d[i] = ca*d[i-1]+cb*v[i-1]+cc*acc[i-1]+cd*acc[i]
        v[i] = cap*d[i-1]+cbp*v[i-1]+ccp*acc[i-1]+cdp*acc[i]
        aa[i] = -ww*d[i]-dcsiw*v[i]

    maxdisp = np.amax(np.absolute(d))
    maxvel = np.amax(np.absolute(v))
    maxacc = np.amax(np.absolute(aa))

    return maxdisp, maxvel, maxacc

def taper(flag, m, samples):
    # m = samples for taper 
    # samples = total samples
    window = kaiser(2*m+1, beta=14)

    if flag == 'front':
        # cut and replace the second half of window with 1s
        ones = np.ones(samples-m-1)
        window = window[0:(m+1)] 
        window = np.concatenate([window, ones]) 

    elif flag == 'end':
        # cut and replace the first half of window with 1s
        ones = np.ones(samples-m-1)
        window = window[(m+1):] 
        window = np.concatenate([ones, window])

    elif flag == 'all':
        ones = np.ones(samples-2*m-1)
        window = np.concatenate([window[0:(m+1)], ones, window[(m+1):]])

    # avoid concatenate error 
    if window.size < samples:
        window = np.append(window, 1)

    if window.size != samples:
        print window.size 
        print samples
        print "[ERROR]: taper and data do not have the same number of samples."
        window = np.ones(samples)
    
    return window 


def seism_appendzeros(flag, t_diff, m, signal):
    """adds zeros in the front and/or at the end of an numpy array
    apply taper before adding
    """
    # if not isinstance(signal, seism_psignal):
    #     return signal

    num = int(t_diff/signal.dt) 
    zeros = np.zeros(num)

    if flag == 'front':
        # applying taper in the front 
        if m != 0: 
            window = taper('front', m, signal.samples)
            # print signal.samples 
            # print signal.accel.size 
            # print window.size 
            signal.accel = signal.accel*window 
            signal.velo = signal.velo*window 
            signal.displ = signal.displ*window 

        # adding zeros in front of data 
        signal.accel = np.append(zeros, signal.accel)
        signal.velo = np.append(zeros, signal.velo)
        signal.displ = np.append(zeros, signal.displ)

    elif flag == 'end':
        if m != 0: 
            # applying taper in the front 
            window = taper('end', m, signal.samples)
            signal.accel = signal.accel*window 
            signal.velo = signal.velo*window 
            signal.displ = signal.displ*window 

        signal.accel = np.append(signal.accel, zeros)
        signal.velo = np.append(signal.velo, zeros)
        signal.displ = np.append(signal.displ, zeros)

    signal.samples += num 
    return signal
# end of seism_appendzeros

def seism_cutting(flag, t_diff, m, signal):
    """cut data in the front or at the end of an numpy array 
    apply taper after cutting 
    """
    # if not isinstance(signal, seism_psignal):
    #     return signal

    num = int(t_diff/signal.dt) 

    if flag == 'front':
        signal.accel = signal.accel[num:]
        signal.velo = signal.velo[num:]
        signal.displ = signal.displ[num:]
        signal.samples -= num 

        # applying taper in the front 
        window = taper('front', m, signal.samples)
        signal.accel = signal.accel*window 
        signal.velo = signal.velo*window 
        signal.displ = signal.displ*window 


    elif flag == 'end': 
        num *= -1 
        # print num 
        # print signal.accel.size 
        signal.accel = signal.accel[:num]
        signal.velo = signal.velo[:num]
        signal.displ = signal.displ[:num]
        signal.samples += num 
        # print signal.samples 
        # print signal.accel.size 

        # applying taper in the front 
        window = taper('end', m, signal.samples)
        signal.accel = signal.accel*window 
        signal.velo = signal.velo*window 

    return signal
# end of seism_cutting

