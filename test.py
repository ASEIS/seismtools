#!/usr/bin/env python
import os
def search_file(dirname, info):
	file_dir = {'HN':'', 'V1':'', 'BH':'', 'V2':''}
	for fp in os.listdir(dirname):
		if (info in fp) and ('HN' in fp): 
			file_dir['HN'] = fp 
			print fp 
			print file_dir.values()


search_file('1', 'CI')
