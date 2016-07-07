#!/usr/bin/env python
"""
# =============================================================================
# The program processes data reference and simulation signals.
# Including rotation, making equal dt, and synchronization
# of starting and ending time.
# =============================================================================
"""
from __future__ import division, print_function
import math
import numpy as np
from scipy import interpolate
from seism import seism_psignal, s_filter
from stools import seism_cutting, seism_appendzeros
# import matplotlib.pyplot as plt

def reverse_up_down(station):
    """
    reverse up down component
    """
    # station has 3 components [ns, ew, ud]
    # only need to flip the 3rd one
    station[2].accel *= -1
    station[2].velo *= -1
    station[2].displ *= -1

    return station
# end of reverse_up_down

def scale_from_m_to_cm(station):
    # scales timeseries from meters to centimeters
    for i in range(0, len(station)):
        station[i].accel *= 100
        station[i].velo *= 100
        station[i].displ *= 100

    return station
# end of scale_data

def get_azimuth():
    """
    Get the azimuth for rotation from user.
    """
    azimuth = ''
    while not azimuth:
        a = raw_input('== Enter azimuth for rotation (optional): ')

        # if user choose not to rotate
        if not a:
            return azimuth

        try:
            azimuth = float(a)
        except ValueError:
            print("[ERROR]: invalid azimuth.")
    return azimuth
# end of get_azimuth

def rotate(station, azimuth):
    """
    Station = [psignal_ns, psignal_ew, psignal_up]
    """
    # checking instance
    if len(station) != 3:
        return station
    for s in station:
        if not isinstance(s, seism_psignal):
            return station

    psignal_ns = station[0]
    psignal_ew = station[1]
    psignal_up = station[2]

    # Nothing to do
    if not azimuth:
        return station

    # rotate data in North and East
    matrix = np.array([(math.cos(math.radians(azimuth)),
                        -math.sin(math.radians(azimuth))),
                       (math.sin(math.radians(azimuth)),
                        math.cos(math.radians(azimuth)))])
    [psignal_ns.accel, psignal_ew.accel] = matrix.dot([psignal_ns.accel,
                                                       psignal_ew.accel])
    [psignal_ns.velo, psignal_ew.velo] = matrix.dot([psignal_ns.velo,
                                                     psignal_ew.velo])
    [psignal_ns.disp, psignal_ew.disp] = matrix.dot([psignal_ns.displ,
                                                     psignal_ew.displ])

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of rotate

# ============================================================================

def get_dt():
    dt = ''
    while not dt:
        d = raw_input("== Enter common dt of two signals: ")
        try:
            dt = float(d)
        except ValueError:
            print("[ERROR]: invalid dt.")
    return dt
# end of get_dt

def get_fmax():
    fmax = ''
    while not fmax:
        f = raw_input("== Enter the maximum frequency for decimation: ")
        try:
            fmax = float(f)
        except ValueError:
            print("[ERROR]: invalid fmax.")
    return fmax
# end of get_fmax

def interp(data, samples, old_dt, new_dt):
    """
    Call interpolate on given data
    """
    old_t = np.arange(0, samples*old_dt, old_dt)
    if old_t.size == samples+1:
        old_t = old_t[:-1]

    f = interpolate.interp1d(old_t, data, 'linear', bounds_error=False)

    new_t = np.arange(0, samples*old_dt, new_dt)
    new_data = f(new_t)

    # eliminate NaN values
    for i in range(1, new_data.size-1):
        if np.isnan(new_data[i]):
            if not np.isnan(new_data[i+1]):
                new_data[i] = (new_data[i-1] + new_data[i+1])/2
            else:
                new_data[i] = new_data[i-1]

    if np.isnan(new_data[-1]):
        new_data[-1] = new_data[-2]
    # using plot to test
    # plt.plot(t,data,'r',new_t,new_data,'b')
    # plt.show()

    return new_data
# end of interpolate

def process_signal_dt(signal, dt, fmax):
    """
    Processes signal with common dt and fmax.
    """
    # call low_pass filter at fmax
    signal.accel = s_filter(signal.accel, signal.dt, type='lowpass',
                            family='butter', fmax=fmax,
                            N=4, rp=0.1, rs=100)
    signal.velo = s_filter(signal.velo, signal.dt, type='lowpass',
                           family='butter', fmax=fmax,
                           N=4, rp=0.1, rs=100)
    signal.displ = s_filter(signal.displ, signal.dt, type='lowpass',
                            family='butter', fmax=fmax,
                            N=4, rp=0.1, rs=100)

    # interpolate
    signal.accel = interp(signal.accel, signal.samples, signal.dt, dt)
    signal.velo = interp(signal.velo, signal.samples, signal.dt, dt)
    signal.displ = interp(signal.displ, signal.samples, signal.dt, dt)

    signal.samples = signal.accel.size
    signal.dt = dt

    return signal
# end of process

def process_dt(station1, station2, dt, fmax):
    """
    Process all signals in two stations to have common dt
    """
    # process signals in stations
    for i in range(0, 3):
        station1[i] = process_signal_dt(station1[i], dt, fmax)
        station2[i] = process_signal_dt(station2[i], dt, fmax)

    return station1, station2
# end of process_dt

# ============================================================================
def get_earthq():
    """
    Get the earthquake start time
    """
    time = raw_input("== Enter the earthquake start time (#:#:#.#): ")
    time = time.split(':')
    if len(time) < 3:
        print("[ERROR]: invalid time format.")
        return get_earthq()
    else:
        for i in range(0, len(time)):
            try:
                time[i] = float(time[i])
            except ValueError:
                print("[ERROR]: invalid time format.")
                return get_earthq()

    # time = [hour, min, sec, frac]
    return time
# end of get_earthq

def get_leading():
    """
    Get the simulation leading time
    """
    lt = ''
    while not lt:
        t = raw_input("== Enter the simulation leading time (sec): ")
        try:
            lt = float(t)
            return lt
        except ValueError:
            print("[ERROR]: invalid leading time.")
    # return lt
# end of get_leading

def synchronize(station1, station2, stamp, eqtimestamp, leading):
    """
    synchronize the stating time and ending time of data arrays in two signals
    signal1 = data signal; signal2 = simulation signal
    """
    if not stamp:
        return station1, station2

    # time in sec = hr*3600 + min*60 + sec + frac*0.1
    start = stamp[0]*3600 + stamp[1]*60 + stamp[2]
    eq_time = eqtimestamp[0]*3600 + eqtimestamp[1]*60 + eqtimestamp[2]
    sim_start = eq_time - leading

    for i in range(0, 3):
        signal1 = station1[i]
        signal2 = station2[i]
        samples1 = signal1.samples
        samples2 = signal2.samples

        dt = signal1.dt # same dt of two signals

        # synchronize the start time
        if start < sim_start:
            # data time < sim time < earthquake time; cutting data array
            signal1 = seism_cutting('front', (sim_start - start),
                                    20, signal1, False)

        elif start > eq_time:
            # sim time < earthquake time < data time; adding zeros in front
            signal1 = seism_appendzeros('front', (start - eq_time),
                                        20, signal1)
            signal2 = seism_cutting('front', (eq_time - sim_start),
                                    20, signal2, False)

        else:
            # sim time < data time < earthquake time; adding zeros
            signal1 = seism_appendzeros('front', (start - sim_start),
                                        20, signal1)

        # synchronize the ending time
        data_time = dt * samples1 # total time of data signal
        end = start + data_time
        sim_time = dt * samples2 # total simulation time
        sim_end = sim_start + sim_time

        if sim_end < end:
            # adding zeros in simulation signal
            # signal2 = seism_appendzeros('end', (end - sim_end), 20, signal2)
            signal1 = seism_cutting('end', (end - sim_end), 20, signal1, False)

        elif end < sim_end:
            # cutting from simulation signal
            signal2 = seism_cutting('end', (sim_end - end), 20, signal2, False)
        else:
            pass

        # scale the data if they have one sample in difference after synchronizing
        if signal1.samples == signal2.samples+1:
            seism_appendzeros('end', signal2.dt, 20, signal2)
        elif signal2.samples == signal1.samples+1:
            seism_appendzeros('end', signal1.dt, 20, signal1)

        station1[i] = signal1
        station2[i] = signal2

    return station1, station2
# end of synchronize
