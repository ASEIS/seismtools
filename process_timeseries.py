#!/usr/bin/env python
"""
# ==============================================================================
# The program is to read input seismograms; process their signals.
# ==============================================================================
"""
from __future__ import division, print_function
import os
import sys
import argparse

from ptools import get_bands, read_file, print_bbp, filter_data, \
    read_stamp, check_data, read_unit_bbp, synchronize_all_stations
from gof_data_sim import get_dt, get_azimuth, get_leading, get_earthq, \
    get_fmax, rotate, scale_from_m_to_cm, process_signal_dt, reverse_up_down

def get_out():
    """
    Get the path of output directory from user
    """
    outdir = ''

    # get the destination saving outputs
    while not outdir:
        outdir = raw_input('== Enter name of the directory to store outputs: ')

    # check existence of target directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return outdir

def process_station_dt(station, common_dt, fmax):
    """
    Process the station to set a common dt
    """
    for i in range(0, 3):
        station[i] = process_signal_dt(station[i],
                                       common_dt,
                                       fmax)
    return station

def process(obs_file, obs_data, stations, params):
    """
    This method processes the signals in each pair of stations.
    Processing consists on scaling, rotation, decimation, alignment
    and other things to make both signals compatible to apply GOF method.
    obs_data: recorded data
    stations: simulation
    """
    # rotate synthetics
    stations = [rotate(station, params['azimuth']) for station in stations]

    # process signals to have the same dt
    if obs_data is not None:
        obs_data = process_station_dt(obs_data,
                                      params['commondt'],
                                      params['decifmax'])
    stations = [process_station_dt(station,
                                   params['commondt'],
                                   params['decifmax']) for station in stations]

    # Read obs_file timestamp if needed
    stamp = None
    if obs_data is not None:
        stamp = read_stamp(obs_file)

    # synchronize starting and ending time of data arrays
    obs_data, stations = synchronize_all_stations(obs_data,
                                                  stations,
                                                  stamp,
                                                  params['eq_time'],
                                                  params['leading'])

    # Check number of samples
    if obs_data is not None:
        num_samples = obs_data[0].samples
    else:
        num_samples = stations[0][0].samples

    for station in stations:
        if station[0].samples != num_samples:
            print("[ERROR]: two timseries do not have the same number"
                  " of samples after processing.")
            sys.exit(-1)

    # Check the data
    if obs_data is not None:
        if not check_data(obs_data):
            print("[ERROR]: processed recorded data contains errors!")
            sys.exit(-1)
    for station in stations:
        if not check_data(station):
            print("[ERROR]: processed simulated data contains errors!")
            sys.exit(-1)

    # All done
    return obs_data, stations
# end of process

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Processes a numer of "
                                     "timeseries files and prepares them "
                                     "for plotting.")
    parser.add_argument("--obs", dest="obs_file",
                        help="input file containing recorded data")
    parser.add_argument("--leading", type=float, dest="leading",
                        help="leading time for the simulation (seconds)")
    parser.add_argument("--eq-time", dest="eq_time",
                        help="earthquake start time (HH:MM:SS.CCC)")
    parser.add_argument("--azimuth", type=float, dest="azimuth",
                        help="azimuth for rotation (degrees)")
    parser.add_argument("--dt", type=float, dest="commondt",
                        help="common dt for the signals")
    parser.add_argument("--decimation-freq", type=float, dest="decifmax",
                        help="maximum frequency for decimation")
    parser.add_argument("--bands", dest="bands",
                        help="sequence of sample rates")
    parser.add_argument("--output-dir", dest="outdir",
                        help="output directory for the outputs")
    parser.add_argument('input_files', nargs='*')
    args = parser.parse_args()

    # Input files
    files = args.input_files
    obs_file = args.obs_file

    if len(files) < 1 or len(files) == 1 and obs_file is None:
        print("[ERROR]: Please provide at least two timeseries to process!")
        sys.exit(-1)

    # Ask user for any missing input parameters
    params = {}
    if args.outdir is None:
        params['outdir'] = get_out()
    else:
        params['outdir'] = args.outdir
    if args.bands is None:
        params['bands'] = get_bands()
    else:
        freqs = args.bands.replace(',', ' ').split()
        if len(freqs) < 2:
            print("[ERROR]: Invalid frequencies!")
            sys.exit(-1)
        try:
            freqs = [float(freq) for freq in freqs]
        except ValueError:
            print("[ERROR]: Invalid frequencies")
        for i in range(0, len(freqs)-1):
            if freqs[i] >= freqs[i+1]:
                print("[ERROR]: Invalid sequence of sample rates")
        params['bands'] = freqs
    if args.decifmax is None:
        params['decifmax'] = get_fmax()
    else:
        params['decifmax'] = args.decifmax
    if args.commondt is None:
        params['commondt'] = get_dt()
    else:
        params['commondt'] = args.commondt
    if args.azimuth is None:
        params['azimuth'] = get_azimuth()
    else:
        params['azimuth'] = args.azimuth
    if args.eq_time is None:
        params['eq_time'] = get_earthq()
    else:
        tokens = args.eq_time.split(':')
        if len(tokens) < 3:
            print("[ERROR]: Invalid time format!")
            sys.exit(-1)
        try:
            params['eq_time'] = [float(token) for token in tokens]
        except ValueError:
            print("[ERROR]: Invalid time format!")
            sys.exit(-1)
    if args.leading is None:
        params['leading'] = get_leading()
    else:
        params['leading'] = args.leading

    return obs_file, files, params

def read_files(obs_file, input_files):
    """
    Reads all input files
    """
    # read obs data
    obs_data = None
    if obs_file is not None:
        obs_data = read_file(obs_file)
        # Make sure we got it
        if not obs_data:
            print("[ERROR]: Reading obs file: %s!" % (obs_file))
            sys.exit(-1)
        # Fix units if needed
        if obs_file.lower().endswith(".bbp"):
            units = read_unit_bbp(obs_file)
            # If in meters, scale to cm
            if units == "m":
                obs_data = scale_from_m_to_cm(obs_data)
        else:
            # In Hercule files, for observation files data is already in
            # cm, so nothing to do here!
            #obs_data = reverse_up_down(obs_data)
            pass

    # reads signals
    stations = []
    for input_file in input_files:
        station = read_file(input_file)
        # Make sure we got it
        if not station:
            print("[ERROR]: Reading input file: %s!" % (input_file))
            sys.exit(-1)
        # Fix units if needed
        if input_file.lower().endswith(".bbp"):
            units = read_unit_bbp(input_file)
            # If in meters, scale to cm
            if units == "m":
                station = scale_from_m_to_cm(station)
        else:
            # Hercule file, need to scale and flip up/down component
            station = scale_from_m_to_cm(station)
            station = reverse_up_down(station)
        stations.append(station)

    # all done
    return obs_data, stations

def process_main():
    """
    Main function for preparing seismograms for plotting
    """
    # First let's get all aruments that we need
    obs_file, input_files, params = parse_arguments()
    obs_data, stations = read_files(obs_file, input_files)

    # processing signals
    obs_data, stations = process(obs_file, obs_data, stations, params)

    # final filtering step
    if obs_data is not None:
        for i in range(0, 3):
            obs_data[i] = filter_data(obs_data[i],
                                      params['bands'][0],
                                      params['bands'][-1])
    for station in stations:
        for i in range(0, 3):
            station[i] = filter_data(station[i],
                                     params['bands'][0],
                                     params['bands'][-1])

    # write processed files
    if obs_data is not None:
        obs_file_out = os.path.join(params['outdir'],
                                    "p-%s" % obs_file.split('/')[-1])
        print_bbp(obs_file, obs_file_out, obs_data)

    for input_file, station in zip(input_files, stations):
        out_file = os.path.join(params['outdir'],
                                "p-%s" % input_file.split('/')[-1])
        print_bbp(input_file, out_file, station)
# end of process_main

# ============================ MAIN ==============================
if __name__ == "__main__":
    process_main()
# end of main program
