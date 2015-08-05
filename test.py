#!/usr/bin/env python
def update_labels(labels):
	# add labels for the parameters used to calculate scores
	o = ['_', '_NS_', '_EW_', '_UD_']
	p = ['PGD', 'PGV', 'PGA', 'A', 'E', 'DUR']
	d = ['D', 'S']
	p_labels = []

	for i in range(0, len(p)):
		for j in range(0, len(d)):
			for k in range(0, len(o)):
				p_labels.append(p[i] + o[k] + d[j])
	print p_labels

	pass 

update_labels([])