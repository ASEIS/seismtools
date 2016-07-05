#!/usr/bin/env python
"""
Utility to convert Hercules time history files to BBP format
"""
from __future__ import division, print_function

# Import python modules
import os
import sys

def write_bbp_header(out_fp, file_type, file_unit,
                     lon, lat, station, time):
    """
    This function writes the bbp header
    """
    # Write header
    out_fp.write("# Station: %s\n" % (station))
    out_fp.write("#    time= %s\n" % (time))
    out_fp.write("#     lon= %s\n" % (lon))
    out_fp.write("#     lat= %s\n" % (lat))
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
    if len(sys.argv) < 3:
        print("Usage: hercules2bbp input_hercules_file output_bbp_file_stem")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file_dis = "%s.dis.bbp" % (sys.argv[2])
    output_file_vel = "%s.vel.bbp" % (sys.argv[2])
    output_file_acc = "%s.acc.bbp" % (sys.argv[2])

    # Covert from Hercules to BBP format by selecting the velocity data
    ifile = open(input_file)
    o_dis_file = open(output_file_dis, 'w')
    o_vel_file = open(output_file_vel, 'w')
    o_acc_file = open(output_file_acc, 'w')
    write_bbp_header(o_dis_file, "displacement", 'm', 0.0, 0.0,
                     "NoName", "0:0:0.0")
    write_bbp_header(o_vel_file, "velocity", 'm/s', 0.0, 0.0,
                     "NoName", "0:0:0.0")
    write_bbp_header(o_acc_file, "acceleration", 'm/s^2', 0.0, 0.0,
                     "NoName", "0:0:0.0")
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
