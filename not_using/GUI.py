#!/usr/bin/env python
# try:
#     fp = open("nother")
# except IOError, e:
#     # print e.errno
#     print e


from Tkinter import * 
import tkMessageBox
import Tkinter as tk 
from main import read_list
import sys
import os
from smc import *
from compare_signals import *

file_list = []
box = []
root = tk.Tk()

def get_filename(): 
	global file_list 
	file_list = []
	global box 
	# if filenam is not given with command 
	if len(sys.argv) == 1:
		flist = E1.get()
		file_list = flist.split()

	# one or more filenames are given; get a list of them
	else: 
		file_list = sys.argv[1:]
	print file_list
	read_list(file_list)
	print file_list

	for f in file_list: 
	    var = tk.IntVar()
	    l = Checkbutton(root, text=f, variable=var)
	    box.append((f, var))
	    l.pack(padx=5, pady=10, side=TOP, anchor = 'w')

def check_state():
	for text, var in box: 
		if var.get():
			print text
			# print "fdfdhfhd"
		else: 
			print "dfdd"



def compare():
   tkMessageBox.showinfo( "Hello Python", "Hello World")


L1 = Label(root, text="file\directory name: ")
L1.pack(side = LEFT)
E1 = Entry(root, bd =5)

E1.pack(side = RIGHT)
add = tk.Button(root, text ="ADD + ", command = get_filename)
add.pack(padx=5, pady=10, side = RIGHT)



B1 = tk.Button(root, text ="Compare .txt", command = check_state)
B2 = tk.Button(root, text ="Compare .her", command = compare)
B3 = tk.Button(root, text ="Generate .txt", command = compare)
B4 = tk.Button(root, text ="Generate .her", command = compare)


B1.pack(padx=5, pady=10, side=LEFT)
B2.pack(padx=5, pady=20, side=LEFT)
B3.pack(padx=5, pady=20, side=LEFT)
B4.pack(padx=5, pady=20, side=LEFT)



root.mainloop()







