#!/usr/bin/env python
# ==========================================
# The program contians the classes holding data.
# ==========================================
import numpy as np
class Data(object):
	def __init__(self, dtype, dataX, dataY, dataZ):
		self.dtype = dtype
		self.dataX = dataX
		self.dataY = dataY
		self.dataZ = dataZ
	# end of init
# end of Data

class PlaneData(object):
	def __init__(self, data, pre_dis, pre_vel):
		self.dis = data
		self.pre_dis = pre_dis
		self.pre_vel = pre_vel
	# end of __init__

	def update(self, data, dtype):
		if dtype == 'd':
			pass
		elif dtype == 'v':
			self.pre_dis = data
		elif dtype == 'a':
			self.pre_vel = data
	# end of update
# end pf planeData