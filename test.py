#!/usr/bin/env python
from numpy import sin, cos, pi, linspace
from numpy.random import randn
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

from matplotlib.pyplot import plot, legend, show, hold, grid, figure, savefig


# Generate a noisy signal to be filtered.
t = linspace(-1, 1, 201)
x = (sin(2 * pi * 0.75 * t*(1-t) + 2.1) + 0.1*sin(2 * pi * 1.25 * t + 1) +
    0.18*cos(2 * pi * 3.85 * t))
xn = x + randn(len(t)) * 0.08

# Create an order 3 lowpass butterworth filter.
b, a = butter(10, 0.05, 'low')
# b, a = butter(3, 0.05)
# scipy.signal.lp2hp(b, a, wo=1.0)[source]
# b, a = signal.lp2hp(b, a, wo=0.05)

# # Apply the filter to xn.  Use lfilter_zi to choose the initial condition
# # of the filter.
# zi = lfilter_zi(b, a)
# z, _ = lfilter(b, a, xn, zi=zi*xn[0])

# # Apply the filter again, to have a result filtered at an order
# # the same as filtfilt.
# z2, _ = lfilter(b, a, z, zi=zi*z[0])

# Use filtfilt to apply the filter.
y = filtfilt(b, a, xn)


# Make the plot.
figure(figsize=(10,5))
hold(True)
plot(t, xn, 'b', linewidth=1.75, alpha=0.75)
# plot(t, z, 'r--', linewidth=1.75)
# plot(t, z2, 'r', linewidth=1.75)
plot(t, y, 'k', linewidth=1.75)
legend(('noisy signal',
        'filtfilt'),
        loc='best')
hold(False)
grid(True)
show()
#savefig('plot.png', dpi=65)
# from scipy.io import loadmat
# from scipy.signal import butter, filtfilt
# from matplotlib.pyplot import plot

# signaldata = loadmat('signaldata.mat')

# input_signal = signaldata['input_signal'][0]

# passband = [0.75*2/30, 5.0*2/30]
# b, a = butter(5, passband, 'bandpass')

# y = filtfilt(b, a, input_signal)
# plot(y)
# print 270 in [180, 90, -180]
# import os.path
# filename = "NC.NHC.E.txt"
# num = 0 
# while os.path.isfile(filename):
# 	num += 1
# 	filename = filename[:-3] + str(num) + filename[-4:]
# 	print filename
# 	num += 1
# 	filename = filename[:-3] + str(num) + filename[-4:]
# 	print filename
# 	num += 1
# 	filename = filename[:-3] + str(num) + filename[-4:]
# 	print filename

# # NC.NHC.E.1.txt
# # NC.NHC.E.1.2.txt
# # NC.NHC.E.1.2.3.txt

# open("testing.txt", 'w').close()
# open("testing.txt", 'w').close()