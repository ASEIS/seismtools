#!/usr/bin/env python
# ==========================================
# The program contians the userInput classes that may be used by other program.
# ==========================================
import os
import sys
from pylab import *
from htools import *
class Input(object):
	def __init__(self, *args):
		self.fp = ""
		self.plotType = ""
		self.deltaT = 0.0
		self.simulationTime, self.alongStrike, self.downDip, self.stepAlongStrike, self.stepDownDip = 0.0, 0, 0, 0, 0
		self.magSelect, self.magnitude = [], False
		self.cumulative = True
		self.snapshots, self.numSnapshots = 's', 0
		self.barChoice, self.barMin, self.barMax = False, 0.0, 0.0
		self.colorMap = 'hot'
		self.printDat = False

		if len(args) == 1:
			args = self.read_parameter(args[0])

		filename =''
		if len(args) == 0:
			filename = self.get_fp()
			self.get_plotType()
			self.get_deltaT()
			# self.get_int('simulationTime')
			self.get_sim_time()
			self.get_int('alongStrike')
			self.get_int('downDip')
			# self.check_size(filename)
			check_size(self, filename, self.alongStrike, self.downDip)

			self.get_int('stepAlongStrike')
			self.get_int('stepDownDip')
			self.get_magSelect()
			if len(self.magSelect) == 1:
				self.get_magnitude()

			self.get_cum()
			self.get_snap()
			if self.snapshots == 'm':
				self.get_numsnap()
			self.get_bar()
			if self.barChoice:
				self.get_min_max()
			self.get_colors()
			self.get_printDat()

		elif len(args) > 0: # if parameters are given with command
			args = list(args)
			try:
				filename = self.set_fp(args[0])
				self.set_plotType(args[1])
				self.set_deltaT(args[2])
				# self.set_int_fields('simulationTime', args[3])
				self.set_sim_time(args[3])
				self.set_int_fields('alongStrike', args[4])
				self.set_int_fields('downDip', args[5])
				# self.check_size(filename)
				check_size(self, filename, self.alongStrike, self.downDip)
				self.set_int_fields('stepAlongStrike', args[6])
				self.set_int_fields('stepDownDip', args[7])
				self.set_mag(args[8])
				if len(self.magSelect) == 1:
					self.set_magnitude(args[9])
					del args[9]


				# self.set_scale(args[9])
				self.set_cum(args[9])
				self.set_snap(args[10])
				if self.snapshots == 'm':
					self.set_numsnap(args[11])
					del args[11]

				self.set_bar(args[11])
				if self.barChoice:
					self.set_min_max(args[12:14])
					del args[12:14]

				self.set_colors(args[12])
				self.set_printDat(args[13])
			except IndexError:
				print "[ERROR]: invalid parameters."
				sys.exit()
	# end of __init__

	def read_parameter(self, filename):
		"""loads the file contains parameters"""
		try:
			fp = open(filename, 'r')
			for line in fp:
				params = line.split(' ')
				return params
		except IOError:
			print "[ERROR]: file not found."
			sys.exit()
	# end of read_parameter

	def set_fp(self, filename):
		try:
			fp = open(filename, 'r')
			self.fp = fp
			return filename
		except IOError:
			print "[ERROR]: file not found."
			return False
	# end of set_fp

	def get_fp(self):
		fp = ''
		while not fp:
			fp = raw_input("== Enter the filename of plane data: ")
			fp = self.set_fp(fp)
		return fp
	# end of get_fp

	def set_plotType(self, ptype):
		type_list = ['a', 'v', 'd']
		try:
			ptype = ptype.lower()[0]
			if ptype not in type_list:
				print "[ERROR]: invalid plot type."
				return False
			else:
				self.plotType = ptype
				return True
		except IndexError:
				print "[ERROR]: invalid plot type."
				return False
	# end of set_plotType

	def get_plotType(self):
		ptype = ''
		while not ptype:
			p = raw_input("== Select displacement, velocity, or acceleration (d/v/a) for plot: ")
			ptype = self.set_plotType(p)
	# end of get_plotType

	def set_deltaT(self, dt):
		try:
			dt = float(dt)
		except ValueError:
			print "[ERROR]: invalid deltaT; float ONLY."
			return False

		if dt <= 0.0:
			print "[ERROR]: invalid deltaT; positive float ONLY."
			return False
		self.deltaT = dt
		return True
	# end of set_deltaT


	def get_deltaT(self):
		dt = 0.0
		while not dt:
			dt = raw_input("== Enter dt: ")
			dt = self.set_deltaT(dt)
	# end of get_deltaT

	def set_int_fields(self, varName, value):
		try:
			value = int(value)
			self.__dict__[varName] = value
			return value

		except ValueError:
			print "[ERROR]: invalid value for " + varName + "."
			return 0
	# end of set_int_fields

	def get_sim_time(self):
		sim = 0.0
		while not sim:
			sim = raw_input("== Enter the total simulation time: ")
			sim = self.set_sim_time(sim)
	# end of get_sim_time

	def set_sim_time(self, sim):
		try:
			sim = float(sim)
			self.simulationTime = sim
			return True
		except ValueError:
			print "[ERROR]: invalid value for simulationTime."
			return False
	# end of set_sim_time

	def get_int(self, varName):
		"""get simulationTime, alongStrike, downDip, stepAlongStrike, and stepDownDip from user."""
		value = 0
		while not value:
			value = raw_input("== Enter the value of " + varName + " : ")
			value = self.set_int_fields(varName, value)
	# end of get_int

	def set_mag(self, mag):
		magSelect = []
		magList = ['x','y','z']
		mag = mag.lower().strip()
		for m in mag:
			if m in magList:
				magSelect.append(m)
			elif m not in [' ', ',']:
				print "[ERROR]: invalid axis; x/y/z or their combinations ONLY."
				return []

		if not magSelect:
			print "[ERROR]: invalid axis; x/y/z or their combinations ONLY."
			return []

		magSelect = sorted(list(set(magSelect))) # sort and remove duplicates
		self.magSelect = magSelect
		return magSelect
	# end of set_mag

	def get_magSelect(self):
		magSelect = []
		while not magSelect:
			magSelect = raw_input("== Enter the axis to plot (x/y/z/xy..): ")
			magSelect = self.set_mag(magSelect)

	def set_magnitude(self, magnitude):
		return self.check_bool('magnitude', magnitude)
	# end of set_magnitude

	def get_magnitude(self):
		# decide whether or not to plot magnitude
		m = ''
		while not m:
			m = raw_input("== Would you like to plot the magnitude? (y/n) ").lower()
			m = self.set_magnitude(m)

	def set_scale(self, scale):
		scale_list = ["lin", "linear", "log", "logarithmic"]
		if scale in scale_list:
			self.scale = scale
			return scale
		else:
			print "[ERROR]: invalid scale type."
			return ''
	# end of set_scale


	def get_scale(self):
		"""get the plotting scale from user"""
		scale = ''
		while not scale:
			scale = raw_input("== Select linear or logarithmic for plot: ").lower()
			scale = self.set_scale(scale)
	# end of get_scale

	def set_cum(self, cum):
		return self.check_bool('cumulative', cum)
	# end of set_cum


	def get_cum(self):
		cum = ''
		while not cum:
			cum = raw_input("== Is this a cumulative plot? (y/n) ").lower()
			cum = self.set_cum(cum)
	# end of get_cum

	def set_snap(self, snapshots):
		if snapshots:
			if snapshots[0] in ['s', 'm']:
				self.snapshots = snapshots[0]
				return snapshots[0]
			else:
				print "[ERROR]: invalid snapshot choice."
				return ''
		else:
			print "[ERROR]: missing snapshot choice."
			return ''
	# end of set_snap

	def get_snap(self):
		snapshots =''
		while not snapshots:
			snapshots = raw_input("== Display a single plot or multiple plots? (s/m) ").lower()
			snapshots = self.set_snap(snapshots)
	# end of get_snap

	def set_numsnap(self, num):
		try:
			numSnapshots = int(num)
			if numSnapshots != 0:
				self.numSnapshots = numSnapshots
				return True
			else:
				print "[ERROR]: numSnapshots cannot be 0. "
				return False
		except ValueError:
			print "[ERROR]: invalid input; int ONLY."
			return False
	# end of set_numsnap

	def get_numsnap(self):
		num = 0
		while not num:
			num = raw_input("== Enter how many seconds in between each plot: ")
			num = self.set_numsnap(num)
	# end of get_numsnap

	def set_bar(self, barChoice):
		return self.check_bool('barChoice', barChoice)
	# end of set_bar

	def get_bar(self):
		barChoice = ''
		while not barChoice:
			barChoice = raw_input("== Set colorbar minimum and maximum (y/n): ").lower()
			barChoice = self.set_bar(barChoice)
	# end of get_bar

	def set_min_max(self, values):
		if len(values) == 2:
			try:
				barMin = float(values[0])
				barMax = float(values[1])
				if barMin < barMax:
					self.barMin = barMin
					self.barMax = barMax
					return True
				else:
					print "[ERROR]: barMax has to be greater than barMin."
					return False
			except IndexError and ValueError:
				print "[ERROR]: enter two floats for min and max values."
				return False
		else:
			print "[ERROR]: enter two floats for min and max values."
			return False
	# end of set_min_max


	def get_min_max(self):
		"""get the minimum and maximum value for colorbar"""
		values = []
		while not values:
			values = raw_input("== Enter the MIN and MAX values for colorbar: ").replace(',', ' ').split()
			values = self.set_min_max(values)
	# end of get_min_max

	def set_colors(self, colorChoice):
		if colorChoice.endswith('.cpt'):
			cdict = gmtColormap(colorChoice)
			cus_cmap = matplotlib.colors.LinearSegmentedColormap('custom_colormap',cdict,256)
			self.colorMap = cus_cmap
			return True
		else:
			# TODO: check existence of color map
			self.colorMap = colorChoice
			return True
	# end of set_colors

	def get_colors(self):
		colorChoice = ''
		while not colorChoice:
			colorChoice = raw_input("== Enter a colormap or the path to custom color file: ")
			colorChoice = self.set_colors(colorChoice)
	# end of get_colors

	def set_printDat(self, printDat):
		return self.check_bool('printDat', printDat)
	# end of set_printDat

	def get_printDat(self):
		printDat = ''
		while not printDat:
			printDat = raw_input("== Would you like the print the matrix data (y/n): ")
			printDat = self.set_printDat(printDat)
	# end of get_printDat

	def check_bool(self, varName, value):
		# check the user inputs associated with boolean values
		bool_dict = {'y':True, 'n':False}
		if value:
			if value[0] in bool_dict.keys():
				self.__dict__[varName] = bool_dict[value[0]]
				return True
			else:
				print "[ERROR]: invalid choice for " + varName + "."
				return False
		else:
			print "[ERROR]: missing choice for " + varName + "."
			return False
	# end of check_bool
# end of input


class DatInput(Input):
	"""subclass of Input; used in datToPlot program."""
	def __init__(self, *args):
		# super(DatInput, self).get_deltaT()
		# super(DatInput, self).get_int('simulationTime')
		if len(args) == 0:
			self.get_path()
			super(DatInput, self).get_int('alongStrike')
			super(DatInput, self).get_int('downDip')
			super(DatInput, self).get_snap()
			if self.snapshots == 'm':
				super(DatInput, self).get_numsnap()
			super(DatInput, self).get_colors()
		else:
			args = list(args)
			try:
				self.set_path(args[0])
				super(DatInput, self).set_int_fields('alongStrike', args[1])
				super(DatInput, self).set_int_fields('downDip', args[2])
				super(DatInput, self).set_snap(args[3])
				if self.snapshots == 'm':
					super(DatInput, self).set_numsnap(args[4])
					del args[4]
				super(DatInput, self).set_colors(args[4])
			except IndexError:
				print "[ERROR]: invalid parameters."
				sys.exit()
	# end of __init__

	def set_path(self, path):
		# TODO: add check
		self.path = path
	# end of set_path

	def get_path(self):
		path = ''
		while not path:
			path = raw_input("== Enter the path to data file(s): ")
			self.set_path(path)
	# end of get_path
# end of DatInput class

class ExtractInput(Input):
	"""subclass of Input; used in extractHer program"""
	def __init__(self, *args):
		if len(args) == 0:
			filename = super(ExtractInput, self).get_fp()
			super(ExtractInput, self).get_deltaT()
			super(ExtractInput, self).get_sim_time()
			super(ExtractInput, self).get_int('alongStrike')
			super(ExtractInput, self).get_int('downDip')
			# super(ExtractInput, self).check_size(filename)
			check_size(self, filename, self.alongStrike, self.downDip)
			super(ExtractInput, self).get_int('stepAlongStrike')
			super(ExtractInput, self).get_int('stepDownDip')
			self.get_coor('x')
			self.get_coor('y')
			# self.get_out()
			self.out_path = get_out()
		else:
			args = list(args)
			try:
				filename = super(ExtractInput, self).set_fp(args[0])
				super(ExtractInput, self).set_deltaT(args[1])
				super(ExtractInput, self).set_sim_time(args[2])
				super(ExtractInput, self).set_int_fields('alongStrike', args[3])
				super(ExtractInput, self).set_int_fields('downDip', args[4])
				# super(ExtractInput, self).check_size(filename)
				check_size(self, filename, self.alongStrike, self.downDip)
				super(ExtractInput, self).set_int_fields('stepAlongStrike', args[5])
				super(ExtractInput, self).set_int_fields('stepDownDip', args[6])
				self.set_coor('x', args[7])
				self.set_coor('y', args[8])
				# self.set_out(args[9])
				self.out_path = set_out(args[9])
			except IndexError:
				print "[ERROR]: invalid parameters."
				sys.exit()
	# end of __init__

	def set_coor(self, flag, coordinate):
		try:
			self.__dict__[flag] = float(coordinate)
			return True
		except ValueError:
			print "[ERROR]: invalid coordinate."
			return False
	# end of set_coor

	def get_coor(self, flag):
		coordinate = False
		while not coordinate:
			coordinate = raw_input("== Enter the " + flag + " coordinate of station: ")
			coordinate = self.set_coor(flag, coordinate)
	# end of get_coor
# end of ExtractInput class

class ResponseInput(Input):
	"""subclass of Input; used in response_spectra_map program"""
	def __init__(self, *args):
		if len(args) == 0:
			filename = super(ResponseInput, self).get_fp()
			super(ResponseInput, self).get_deltaT()
			super(ResponseInput, self).get_sim_time()
			super(ResponseInput, self).get_int('dimensionX')
			super(ResponseInput, self).get_int('dimensionY')
			super(ResponseInput, self).get_int('spaceX')
			super(ResponseInput, self).get_int('spaceY')
			check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
			self.get_component()
			super(ResponseInput, self).get_plotType()
			self.get_period()
			self.get_out_space('X')
			self.get_out_space('Y')
			super(ResponseInput, self).get_bar()
			if self.barChoice:
				super(ResponseInput, self).get_min_max()
			super(ResponseInput, self).get_colors()
			super(ResponseInput, self).get_printDat()
			if self.printDat:
				self.out_path = get_out()

		else:
			args = list(args)
			try:
				filename = super(ResponseInput, self).set_fp(args[0])
				super(ResponseInput, self).set_deltaT(args[1])
				super(ResponseInput, self).set_sim_time(args[2])
				super(ResponseInput, self).set_int_fields('dimensionX', args[3])
				super(ResponseInput, self).set_int_fields('dimensionY', args[4])
				# super(ResponseInput, self).check_size(filename)
				super(ResponseInput, self).set_int_fields('spaceX', args[5])
				super(ResponseInput, self).set_int_fields('spaceY', args[6])
				check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
				self.set_component(args[7])
				super(ResponseInput, self).set_plotType(args[8])
				self.set_period(args[9])
				self.set_out_space(args[10], 'X')
				self.set_out_space(args[11], 'Y')
				super(ResponseInput, self).set_bar(args[12])
				if self.barChoice:
					super(ResponseInput, self).set_min_max(args[13:15])
					del args[13:15]
				super(ResponseInput, self).set_colors(args[13])
				super(ResponseInput, self).set_printDat(args[14])
				if self.printDat:
					self.out_path = set_out(args[15])
			except IndexError:
				print "[ERROR]: invalid number of parameters."
				sys.exit()
	# end of __init__

	def set_component(self, component):
		c = ['x', 'y', 'z']
		try:
			component = component.lower()[0]
		except (AttributeError, IndexError) as e:
			print "[Error]: invalid component selection."
			return False

		if component in c:
			self.component = component
			return True
		else:
			print "[Error]: invalid component selection."
			return False
	# end of set_component

	def get_component(self):
		component = ''
		while not component:
			component = raw_input("== Enter the component choice [x/y/z]: ")
			component = self.set_component(component)
	# end of get_component

	def set_period(self, period):
		try:
			self.period = float(period)
			return True
		except ValueError:
			print "[ERROR]: invalid value of period."
			return False
	# end of set_period

	def get_period(self):
		period = 0.0
		while not period:
			period = self.set_period(raw_input("== Enter the value of period: "))
	# end of get_period

	def set_out_space(self, space, axis):
		varName = "outspace" + axis
		try:
			space = int(space)
		except ValueError:
			print "[ERROR]: enter integer value for output space. "

		if axis == 'x':
			sapceChoice = self.spaceX
			dimenChoice = self.dimensionX
		else:
			sapceChoice = self.spaceY
			dimenChoice = self.dimensionY

		# check output space choice with original space and dimension length
		if space%sapceChoice == 0 and dimenChoice%space == 0:
			self.__dict__[varName] = space
			return True
		else:
			print "[ERROR]: output space has to be multiple of original space and be divisible by dimension length."
			return False
	# end of set_out_space

	def get_out_space(self, axis):
		space = 0
		while not space:
			space =  raw_input("== Enter the value of output space in " + axis + " axis: ")
			space = self.set_out_space(space, axis)
	# end of get_out_space
# end of ResponseInput class

class ResponseInput2(ResponseInput):
	"""subclass of ResponseInput; used in response_spectra_map2 program"""
	def __init__(self, *args):
		if len(args) == 0:
			filename = super(ResponseInput2, self).get_fp()
			super(ResponseInput2, self).get_deltaT()
			super(ResponseInput2, self).get_sim_time()
			super(ResponseInput2, self).get_int('dimensionX')
			super(ResponseInput2, self).get_int('dimensionY')
			# super(ResponseInput2, self).check_size(filename)
			super(ResponseInput2, self).get_int('spaceX')
			super(ResponseInput2, self).get_int('spaceY')
			check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
			super(ResponseInput2, self).get_component()
			super(ResponseInput2, self).get_plotType()
			super(ResponseInput2, self).get_period()
			self.get_size()
			super(ResponseInput2, self).get_bar()
			if self.barChoice:
				super(ResponseInput2, self).get_min_max()
			super(ResponseInput2, self).get_colors()
			super(ResponseInput2, self).get_printDat()
			if self.printDat:
				self.out_path = get_out()

		else:
			args = list(args)
			try:
				filename = super(ResponseInput2, self).set_fp(args[0])
				super(ResponseInput2, self).set_deltaT(args[1])
				super(ResponseInput2, self).set_sim_time(args[2])
				super(ResponseInput2, self).set_int_fields('dimensionX', args[3])
				super(ResponseInput2, self).set_int_fields('dimensionY', args[4])
				# super(ResponseInput2, self).check_size(filename)
				super(ResponseInput2, self).set_int_fields('spaceX', args[5])
				super(ResponseInput2, self).set_int_fields('spaceY', args[6])
				check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
				super(ResponseInput2, self).set_component(args[7])
				super(ResponseInput2, self).set_plotType(args[8])
				super(ResponseInput2, self).set_period(args[9])
				self.set_size(args[10])
				super(ResponseInput2, self).set_bar(args[11])
				if self.barChoice:
					super(ResponseInput2, self).set_min_max(args[12:14])
					del args[12:14]
				super(ResponseInput2, self).set_colors(args[12])
				super(ResponseInput2, self).set_printDat(args[13])
				if self.printDat:
					self.out_path = set_out(args[14])
			except IndexError:
				print "[ERROR]: invalid number of parameters."
				sys.exit()
	# end of __init__

	def set_size(self, size):
		try:
			size = int(size)
		except ValueError:
			print "[ERROR]: enter integer number of points to load."
			return False

		if 0 < size <= (self.dimensionX/self.spaceX+1):
			self.size = size
			return True
		else:
			print "[ERROR]: size has to be no greater than number of points in X axis."
			return False
	# end of set_size

	def get_size(self):
		size = 0
		while not size:
			size = raw_input("== Enter number of points to load at once: ")
			size = self.set_size(size)
	# end of get_size
# end of ResponseInput2 class

class FourierInput(ResponseInput2):
	"""subclass of ResponseInput2; used in fourier_spectra_map program"""
	def __init__(self, *args):
		if len(args) == 0:
			filename = super(FourierInput, self).get_fp()
			super(FourierInput, self).get_deltaT()
			super(FourierInput, self).get_sim_time()
			super(FourierInput, self).get_int('dimensionX')
			super(FourierInput, self).get_int('dimensionY')
			# super(FourierInput, self).check_size(filename)
			super(FourierInput, self).get_int('spaceX')
			super(FourierInput, self).get_int('spaceY')
			check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
			super(FourierInput, self).get_component()
			super(FourierInput, self).get_plotType()
			self.get_freq()
			super(FourierInput, self).get_size()
			super(FourierInput, self).get_bar()
			if self.barChoice:
				super(FourierInput, self).get_min_max()
			super(FourierInput, self).get_colors()
			super(FourierInput, self).get_printDat()
			if self.printDat:
				self.out_path = get_out()

		else:
			args = list(args)
			try:
				filename = super(FourierInput, self).set_fp(args[0])
				super(FourierInput, self).set_deltaT(args[1])
				super(FourierInput, self).set_sim_time(args[2])
				super(FourierInput, self).set_int_fields('dimensionX', args[3])
				super(FourierInput, self).set_int_fields('dimensionY', args[4])
				# super(FourierInput, self).check_size(filename)
				super(FourierInput, self).set_int_fields('spaceX', args[5])
				super(FourierInput, self).set_int_fields('spaceY', args[6])
				check_size(self, filename, int(self.dimensionY/self.spaceY)+1, int(self.dimensionX/self.spaceX)+1)
				super(FourierInput, self).set_component(args[7])
				super(FourierInput, self).set_plotType(args[8])
				self.set_freq(args[9])
				super(FourierInput, self).set_size(args[10])
				super(FourierInput, self).set_bar(args[11])
				if self.barChoice:
					super(FourierInput, self).set_min_max(args[12:14])
					del args[12:14]
				super(FourierInput, self).set_colors(args[12])
				super(FourierInput, self).set_printDat(args[13])
				if self.printDat:
					self.out_path = set_out(args[14])
			except IndexError:
				print "[ERROR]: invalid number of parameters."
				sys.exit()
	# end of __init__
	def set_freq(self, frequence):
		try:
			self.frequence = float(frequence)
			return True
		except ValueError:
			print "[ERROR]: invalid value of frequence."
			return False
	# end of set_freq

	def get_freq(self):
		frequence = 0.0
		while not frequence:
			frequence = self.set_freq(raw_input("== Enter the value of frequence: "))
	# end of get_freq
# end of FourierInput



def set_out(out_path):
	if '/' in out_path:
		tmp = out_path.split('/')
		out_dir = '/'.join(tmp[:-1])
		if os.path.exists(out_dir):
			return out_path
		else:
			print "[ERROR]: invalid path; directory does not exist."
			return False
	else:
		current = os.getcwd()
		return current + '/' + out_path
# end of set_out

def get_out():
	out_path = ''
	while not out_path:
		out_path = raw_input("== Enter the complete path to output file: ")
		out_path = set_out(out_path)
	return out_path
# end of get_out

def check_size(UserInput, filename, alongStrike, downDip):
	"""compare the size of given file and the estimated size"""
	stat = os.stat(filename)
	act_size = stat.st_size
	est_size = alongStrike*downDip*UserInput.simulationTime/UserInput.deltaT*64*3/8

	if act_size != est_size:
		print "[ERROR]: the size of given file doesn't match given parameters."
		sys.exit()
# end of check_size

# ==================TESTING===================
if __name__ == "__main__":
	if len(sys.argv) > 1:
		argument = tuple(sys.argv[1:])
		i = ResponseInput2(*argument)
	else:
		i = ResponseInput2()

	# print i.__dict__
	# d = DatInput()
	# e = ExtractInput()
	print i.__dict__

