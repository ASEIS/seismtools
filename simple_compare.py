#!/usr/bin/env python
"""
# ==============================================================================
# The program is to plot a simple comparison among timeseries files
# ==============================================================================
"""
from __future__ import division, print_function
import os
import sys

from ptools import read_file
from compare_signals import simple_plot, set_parameter

def simple_compare_main():
    """
    Main function for simple_compare
    """
    num_lines = 0
    output_file = ''

    # Parse command line
    if len(sys.argv) == 1:
        while True:
            num_lines = raw_input("==> How many files to plot: ")
            try:
                num_lines = int(num_lines)
                if num_lines > 0:
                    break
            except ValueError:
                pass
    else:
        # First argument is number of timeseries
        num_lines = int(sys.argv[1])

    if len(sys.argv) < 3:
        output_file = raw_input("==> Output filename: ")
    else:
        # Second argument is output filename
        output_file = sys.argv[2]
    # Allow user to specify screen output by using -
    if output_file == "-":
        output_file = ""

    filenames = []
    # Get input files
    if len(sys.argv) >= num_lines + 3:
        # User provided all input file names
        filenames = sys.argv[3:(num_lines + 3)]
        para = sys.argv[(num_lines + 3):]
    else:
        # Need to ask user
        para = sys.argv[3:]
        remaining = num_lines
        while remaining:
            input_file = raw_input("==> Enter input file %d: " %
                                   (num_lines - remaining + 1))
            filenames.append(input_file)
            remaining = remaining - 1

    # Set all other parameters
    parameter = set_parameter(para)

    # Read data
    stations = [read_file(filename) for filename in filenames]
    filenames = [os.path.basename(filename) for filename in filenames]

    # Create plot
    simple_plot(parameter, filenames, stations, output_file=output_file)

# ============================ MAIN ==============================
if __name__ == "__main__":
    simple_compare_main()
# end of main program
