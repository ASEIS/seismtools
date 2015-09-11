#!/usr/bin/env python
import os
# def search_file(dirname, info):
# 	file_dir = {'HN':'', 'V1':'', 'BH':'', 'V2':''}
# 	for fp in os.listdir(dirname):
# 		if (info in fp) and ('HN' in fp):
# 			file_dir['HN'] = fp
# 			print fp
# 			print file_dir.values()


# search_file('1', 'CI')

def set_labels(bands):
	# generate labels for scores file
	o = ['A', 'N', 'E', 'U']
	b = ['CA', 'SA']
	s = ['T', 'A', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']
	b_label = "BB"

	for i in range(1, len(bands)):
		b_label += ',B'+str(i)
	b_label = b_label.split(',')
	b += b_label

	print o
	print b
	print s
	labels = ['#SIGNAL1', 'SIGNAL2', 'X_COOR', 'Y_COOR', 'EPI_DIS']
	for i in range(0, len(o)):
		for k in range(0, len(s)):
			for j in range(0, len(b)):
				labels.append(o[i]+'_'+b[j]+'_'+s[k])
	return labels
# end of set_labels

labels = set_labels([0.1, 0.25, 0.5])
for l in labels:
	print l