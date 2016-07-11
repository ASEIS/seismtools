#!/usr/bin/env python
"""
# =============================================================================
# The program is to call comparison on TWO 10-column files or TWO 1-column files
# =============================================================================
"""
from __future__ import print_function
import sys
from compare_signals import compare_her, compare_txt, set_parameter

def check_type(filename):
    """
	checks the type of file being .her or .txt
	"""
    try:
        with open(filename) as input_file:
            content = input_file.readlines()
            num_col = len(content[2].split())
            if num_col == 1:
                return 'TXT'
            elif num_col == 10:
                return 'HER'
            else:
                print("[ERROR]: Cannot recognize file type.")
                return False

    except IOError as e:
        print(e)
        return False
# end of check_type

def get_filename():
    file1 = ''
    file2 = ''
    para = []

    # get paths of two files
    if len(sys.argv) >= 2:
        file1 = sys.argv[1]

    if len(sys.argv) >= 3:
        file2 = sys.argv[2]

    if len(sys.argv) > 3:
        para = sys.argv[3:]
    parameter = set_parameter(para)

    while not file1:
        file1 = raw_input('== Enter the path of file1: ')

    while not file2:
        file2 = raw_input('== Enter the path of file2: ')
    return parameter, file1, file2
# end of get_filename

if __name__ == "__main__":
    PARAMETER, FILE1, FILE2 = get_filename()

    # checking file types
    if check_type(FILE1) == check_type(FILE2) == 'TXT':
        compare_txt(PARAMETER, FILE1, FILE2)
    elif check_type(FILE1) == check_type(FILE2) == 'HER':
        compare_her(PARAMETER, FILE1, FILE2, False)
