#!/usr/bin/env python
"""
# ==============================================================================
# The program is to read two timeseries; process their signals;
# calculate their scores with different sample rates;
# and generate 3D matrix for scores.
# ==============================================================================
"""
from __future__ import division, print_function

import os
import sys
import copy
import glob
import math
import argparse
import numpy as np

from process_timeseries import read_files, process
from ptools import read_filelist, get_bands, read_stamp, check_data, read_file, \
    print_her
from gof_engine import print_scores, set_labels, set_mlabels, scores_matrix, \
    filter_data, print_matrix, parameter_to_list
from gof_data_sim import get_dt, get_azimuth, get_leading, get_earthq, \
    get_fmax, rotate, synchronize, scale_from_m_to_cm, process_dt, \
    reverse_up_down

np.seterr(divide='ignore', invalid='ignore')

def get_epicenter():
    """
    Get x and y coordinates of epicenter
    """
    epi = ''
    x = 0.0
    y = 0.0
    while not epi:
        epi = raw_input('== Enter the X and Y coordinates of epicenter: ')
        epi = epi.replace(',', ' ')
        epi = epi.split()
        if len(epi) == 2:
            try:
                x = float(epi[0])
                y = float(epi[1])
                return x, y
            except ValueError:
                print("[ERROR]: invalid coordinates.")
                epi = ''

        print("[ERROR]: invalid coordinates.")
        epi = ''
# end of get_epicenter

def get_in():
    """
    Get the path of input directories
    """
    while True:
        indir1 = ''
        indir2 = ''

        while not indir1:
            indir1 = raw_input('== Enter name of 1st input directory: ')

        while not indir2:
            indir2 = raw_input('== Enter name of 2nd input directory: ')

        # check the existence of two directories
        if (not os.path.exists(indir1)) or (not os.path.exists(indir2)):
            print("[ERROR]: input directory does not exist.")
            # Continue in the while True loop...
        else:
            break

    # Got two valid input directories
    return indir1, indir2

def get_out():
    """
    Get the path of output directory and output file from user
    """
    outdir = ''
    outname1 = ''
    outname2 = ''

    # get the destination saving outputs
    while not outdir:
        outdir = raw_input('== Enter name of the directory to store outputs: ')

    # check existence of target directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    while not outname1:
        outname1 = raw_input('== Enter name of scores file: ')

    while not outname2:
        outname2 = raw_input('== Enter name of metrics file: ')

    path1 = os.path.join(outdir, outname1)
    path2 = os.path.join(outdir, outname2)
    return outdir, path1, path2

def search_file(dirname, info):
    """
    Search for files contains given station code and name
    """
    info = info.replace('.', '')
    for fp in os.listdir(dirname):
        tmp = fp.replace('.', '')
        if (info in tmp) and ('HN' in tmp):
            return fp
            # file_dir['HN'] = fp
        elif (info in tmp) and ('V1' in tmp):
            # file_dir['V1'] = fp
            return fp
        elif (info in tmp) and ('BH' in tmp):
            # file_dir['BH'] = fp
            return fp
        elif (info in tmp) and ('V2' in tmp):
            # file_dir['V2'] = fp
            return fp
        # endif case
    # end for

    # was not found, return same info
    return info
# end of search_file

def find_station(input_directory, station_name):
    """
    Looks into input_directory for data belonging to a station whose name
    matches station_name
    """
    filelist = glob.glob("%s/%s[-,_,.]*.vel.bbp" % (input_directory, station_name))
    # Check if we found anything
    if not filelist:
        return None
    # Only one option
    if len(filelist) == 1:
        return filelist[0]
    # Multiple options, need to figure out which one we want
    print("Multiple files to choose from!")
    print(filelist)
    sys.exit(0)

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Calculates the GOF "
                                     "from a pair or timeseries files or "
                                     "from a list of timeseries.")
    parser.add_argument("--obs", dest="obs_file",
                        help="input file containing recorded data")
    parser.add_argument("--syn", action="append", dest="syn_files",
                        help="input files containing synthetic data")
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
    parser.add_argument("--scores", dest="s_path",
                        help="scores file")
    parser.add_argument("--metrics", dest="m_path",
                        help="metrics file")
    parser.add_argument("--list", dest="filelist",
                        help="list containing files to use")
    parser.add_argument("--input_dir", action="append", dest="indirs",
                        help="input directories for the files in the list")
    parser.add_argument("--epicenter_x", type=float, dest="epicenter_x",
                        help="epicenter coordinates")
    parser.add_argument("--epicenter_y", type=float, dest="epicenter_y",
                        help="epicenter coordinates")
    args = parser.parse_args()

    # Parameters from the user
    params = {}

    # Figure out if we have a pair or timeseries or a list of files to use
    if args.filelist is not None and (args.obs_file is not None or
                                      args.syn_files is not None):
        print("[ERROR]: Specify either list or pair or files!")
        sys.exit(-1)

    if args.filelist is not None:
        # List of files is provided
        params["filelist"] = args.filelist
        if args.indirs is None:
            params["indir1"], params["indir2"] = get_in()
        else:
            if len(args.indirs) == 2:
                params["indir1"] = args.indirs[0]
                params["indir2"] = args.indirs[1]
            else:
                print("[ERROR]: Please specify 2 input directories!")
                sys.exit(-1)
    else:
        # A pair of files
        if args.obs_file is not None:
            len_obs = 1
        else:
            len_obs = 0
        if args.syn_files is not None:
            len_syn = len(args.syn_files)
        else:
            len_syn = 0
        # we need two (and only two!) input files
        if len_obs + len_syn != 2:
            print("[ERROR]: Please specify 2 input files with obs_files "
                  "or syn_files!")
            sys.exit(-1)
        params["obs_file"] = args.obs_file
        params["syn_files"] = args.syn_files

    # Ask user for any missing input parameters
    if args.outdir is None or args.s_path is None or args.m_path is None:
        params['outdir'], params['s_path'], params['m_path'] = get_out()
    else:
        params['outdir'] = args.outdir
        params['s_path'] = os.path.join(args.outdir, args.s_path)
        params['m_path'] = os.path.join(args.outdir, args.m_path)
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
    # Optional
    params['epi_x'] = args.epicenter_x
    params['epi_y'] = args.epicenter_y

    return params
#end parse_arguments

def main_gof():
    """
    Main function for GOF code
    """
    # First let's parse all the arguments that we need
    params = parse_arguments()

    if not "filelist" in params:
        # Two file option!
        obs_data, stations = read_files(params['obs_file'],
                                        params['syn_files'])
        # processing signals
        obs_data, stations = process(params['obs_file'], obs_data,
                                     stations, params)

        # Figure out stations 1 & 2
        if obs_data is None:
            station1 = stations[0]
            station2 = stations[1]
        else:
            station1 = obs_data
            station2 = stations[0]

        # Calculate scores matrix
        parameter, matrix, flag = scores_matrix(station1, station2,
                                                params['bands'])

        # Exit if GOF failed
        if not flag:
            print("[ERROR]: Files not processed!")
            sys.exit(-1)

        print_matrix(params['s_path'], matrix)
        # End: Two-Files Option
    else:
        # Start: List of Files Option
        station_list, coor_x, coor_y = read_filelist(params['filelist'])

        # Get coordinates
        if coor_x[0] and coor_y[0]:
            if params["epi_x"] is not None and params["epi_y"] is not None:
                epi_x = params["epi_x"]
                epi_y = params["epi_y"]
            else:
                # Ask user
                epi_x, epi_y = get_epicenter()
        else:
            epi_x = epi_y = 0.0

        # open output files
        try:
            f = open(params['s_path'], 'w')
            m = open(params['m_path'], 'w')
            u = open(os.path.join(params['outdir'], "unprocessed.txt"), 'w')
        except IOError as err:
            print(err)

        # prepares formats for output files
        labels = set_labels(params['bands'])
        m_labels = set_mlabels()

        d = '{:>12}'*2 + '{:>12.8}'*(len(labels)-2) + '\n'
        f.write(d.format(*labels))

        d = '{:>12}'*2 + '{:>12.6}'*(len(m_labels)-2) + '\n'
        m.write(d.format(*m_labels))
        f.close()
        m.close()

        # create empty list for pairs that fail to be processed
        unprocessed = []

        # loop of the list of pairs given in the list-file
        for i in range(0, len(station_list)):
            # capture full path to files
            file1 = find_station(params['indir1'], station_list[i])
            file2 = find_station(params['indir2'], station_list[i])

            if file1 is None or file2 is None:
                # Add to list of unprocessed stations
                if file1 is None:
                    tmpmsg = "%s (no data)" % (station_list[i])
                    unprocessed.append(tmpmsg)
                if file2 is None:
                    tmpmsg = "%s (no synthetic)" % (station_list[i])
                    unprocessed.append(tmpmsg)

                print("...Ignoring station:   %s" % (station_list[i]))
                continue

            # Both files are available, attempts to process...
            print("\n...Processing pair: " + file1 + " - " + file2)

            # computes epicentral distance
            x = coor_x[i]
            y = coor_y[i]
            epdist = math.sqrt((x-epi_x)**2+(y-epi_y)**2)
            coord = [x, y, epdist]

            # reads signals
            obs_data, stations = read_files(file1, [file2])

            # processing signals
            obs_data, stations = process(file1, obs_data,
                                         stations, params)
            station1 = obs_data
            station2 = stations[0]

            if station1 and station2:
                parameter, matrix, flag = scores_matrix(station1,
                                                        station2,
                                                        params['bands'])
                # Sanity check to avoid division by zero
                if not flag:
                    tmpmsg = "%s (div by zero)" % (station_list[i])
                    unprocessed.append(tmpmsg)
                    print("...Ignoring station:   %s (div by zero)" %
                          (station_list[i]))
                    continue
                # end if: sanity check

                parameter = parameter_to_list(parameter)

                # print scores
                print_scores([file1, file2], coord,
                             params['s_path'], [], matrix)
                # print values used to calculate scores
                print_scores([file1, file2], coord, params['m_path'],
                             parameter, np.array([]))
            else:
                pass
            # end if station1 and station2
        # end loop of the list of pairs

        for pair in unprocessed:
            u.write("%s\n" % pair)

    #end of if instance switch

    print("\n[DONE]")

# ============================ MAIN ==============================
if __name__ == "__main__":
    main_gof()
# end of main program
