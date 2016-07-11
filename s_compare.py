#!/usr/bin/env python
"""
# ============================================================================
# The program is a simplied version of compare.py. It plots only velocity for
# data and FAS, and acceleration for Response.
# ============================================================================
"""
from compare import check_type, get_filename
from compare_signals import compare_her

if __name__ == "__main__":
    parameter, file1, file2 = get_filename()
    if not (check_type(file1) == check_type(file2) == 'HER'):
        print "[ERROR]: unsupported file type."
        # return

    compare_her(parameter, file1, file2, True)
