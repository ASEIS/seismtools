#!/usr/bin/env python
"""
Utility to convert Hercules time history files to BBP format
"""
from __future__ import division, print_function

# Import python modules
import os
import argparse

def write_bbp_header(out_fp, file_type, file_unit, args):
    """
    This function writes the bbp header
    """
    # Write header
    out_fp.write("# Station: %s\n" % (args.station_name))
    out_fp.write("#    time= %s\n" % (args.time))
    out_fp.write("#     lon= %s\n" % (args.longitude))
    out_fp.write("#     lat= %s\n" % (args.latitude))
    out_fp.write("#   units= %s\n" % (file_unit))
    out_fp.write("#\n")
    out_fp.write("# Data fields are TAB-separated\n")
    out_fp.write("# Column 1: Time (s)\n")
    out_fp.write("# Column 2: N/S component ground "
                 "%s (%s) (+ is 000)\n" % (file_type, file_unit))
    out_fp.write("# Column 3: E/W component ground "
                 "%s (%s) (+ is 090)\n" % (file_type, file_unit))
    out_fp.write("# Column 4: U/D component ground "
                 "%s (%s) (+ is upward)\n" % (file_type, file_unit))
    out_fp.write("#\n")

def hercules2bbp_main():
    """
    Main function for hercules to bbp converter
    """
    parser = argparse.ArgumentParser(description="Converts a Hercules "
                                     "file to BBP format, generating "
                                     "displacement, velocity and acceleration "
                                     "BBP files.")
    parser.add_argument("-s", "--station-name", dest="station_name",
                        default="NoName",
                        help="provides the name for this station")
    parser.add_argument("--lat", dest="latitude", type=float, default=0.0,
                        help="provides the latitude for the station")
    parser.add_argument("--lon", dest="longitude", type=float, default=0.0,
                        help="provides the longitude for the station")
    parser.add_argument("-t", "--time", default="00/00/00,0:0:0.0 UTC",
                        help="provides timing information for this timeseries")
    parser.add_argument("input_file", help="Hercules input timeseries")
    parser.add_argument("output_stem",
                        help="output BBP filename stem without the "
                        " .{dis,vel,acc}.bbp extensions")
    parser.add_argument("-d", dest="output_dir", default="",
                        help="output directory for the BBP file")
    args = parser.parse_args()

    input_file = args.input_file
    output_file_dis = "%s.dis.bbp" % (os.path.join(args.output_dir,
                                                   args.output_stem))
    output_file_vel = "%s.vel.bbp" % (os.path.join(args.output_dir,
                                                   args.output_stem))
    output_file_acc = "%s.acc.bbp" % (os.path.join(args.output_dir,
                                                   args.output_stem))

    # Covert from Hercules to BBP format by selecting the velocity data
    ifile = open(input_file)
    o_dis_file = open(output_file_dis, 'w')
    o_vel_file = open(output_file_vel, 'w')
    o_acc_file = open(output_file_acc, 'w')
    write_bbp_header(o_dis_file, "displacement", 'm', args)
    write_bbp_header(o_vel_file, "velocity", 'm/s', args)
    write_bbp_header(o_acc_file, "acceleration", 'm/s^2', args)
    for line in ifile:
        line = line.strip()
        # Skip comments
        if line.startswith("#") or line.startswith("%"):
            pieces = line.split()[1:]
            # Write header
            if len(pieces) >= 10:
                o_dis_file.write("# her header: # %s %s %s %s\n" %
                                 (pieces[0], pieces[1], pieces[2], pieces[3]))
                o_vel_file.write("# her header: # %s %s %s %s\n" %
                                 (pieces[0], pieces[4], pieces[5], pieces[6]))
                o_acc_file.write("# her header: # %s %s %s %s\n" %
                                 (pieces[0], pieces[7], pieces[8], pieces[9]))
            else:
                o_dis_file.write("# her header: %s\n" % (line))
            continue
        pieces = line.split()
        pieces = [float(piece) for piece in pieces]
        # Write timeseries to files. Please not that Hercules files have
        # the vertical component positive pointing down so we have to flip it
        # here to match the BBP format in which vertical component points up
        o_dis_file.write("%1.9E %1.9E %1.9E %1.9E\n" %
                         (pieces[0], pieces[1], pieces[2], -1 * pieces[3]))
        o_vel_file.write("%1.9E %1.9E %1.9E %1.9E\n" %
                         (pieces[0], pieces[4], pieces[5], -1 * pieces[6]))
        o_acc_file.write("%1.9E %1.9E %1.9E %1.9E\n" %
                         (pieces[0], pieces[7], pieces[8], -1 * pieces[9]))

    # All done, close everything
    ifile.close()
    o_dis_file.close()
    o_vel_file.close()
    o_acc_file.close()

# ============================ MAIN ==============================
if __name__ == "__main__":
    hercules2bbp_main()
# end of main program
