#!/usr/bin/env python
# ===========================================Example One============================================
# from numpy import sin, cos, pi, linspace
# from numpy.random import randn
# from scipy.signal import lfilter, lfilter_zi, filtfilt, butter, ellip

# from matplotlib.pyplot import plot, legend, show, hold, grid, figure, savefig


# # Generate a noisy signal to be filtered.
# t = linspace(-1, 1, 201)
# x = (sin(2 * pi * 0.75 * t*(1-t) + 2.1) + 0.1*sin(2 * pi * 1.25 * t + 1) +
#     0.18*cos(2 * pi * 3.85 * t))
# xn = x + randn(len(t)) * 0.08

# # Create an order 3 lowpass butterworth filter.
# # b, a = butter(3, 0.05)
# # b, a = butter(10, 0.05, 'low')

#  # scipy.signal.ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba')
# # b, a = ellip(1, 1, 1, 1, 'high', analog=True)
# b, a = ellip(N = 17, rp = 0.01, rs = 60, Wn = 0.05/((1.0/0.05)/2.0), btype = 'highpass', analog=False)


# # # Apply the filter to xn.  Use lfilter_zi to choose the initial condition
# # # of the filter.
# # zi = lfilter_zi(b, a)
# # z, _ = lfilter(b, a, xn, zi=zi*xn[0])

# # # Apply the filter again, to have a result filtered at an order
# # # the same as filtfilt.
# # z2, _ = lfilter(b, a, z, zi=zi*z[0])

# # Use filtfilt to apply the filter.
# y = filtfilt(b, a, xn)


# # Make the plot.
# figure(figsize=(10,5))
# hold(True)
# plot(t, xn, 'b', linewidth=1.75, alpha=0.75)
# # plot(t, z, 'r--', linewidth=1.75)
# # plot(t, z2, 'r', linewidth=1.75)
# plot(t, y, 'k', linewidth=1.75)
# legend(('noisy signal',
#         'filtfilt'),
#         loc='best')
# hold(False)
# grid(True)
# show()

# ===========================================Example Two============================================
import pylab
import scipy
import scipy.signal
# [b,a] = scipy.signal.ellip(6,3,50,300.0/500.0);
[b,a] = scipy.signal.ellip(N = 17, rp = 0.01, rs = 60, Wn = 0.05/((1.0/0.05)/2.0), btype = 'highpass', analog=False);

import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
plt.title('Digital filter frequency response')
ax1 = fig.add_subplot(111)
h,w = scipy.signal.freqz(b, a)
plt.semilogy(h, np.abs(w), 'b')
plt.semilogy(h, abs(w), 'b')
plt.ylabel('Amplitude (dB)', color='b')
plt.xlabel('Frequency (rad/sample)')
plt.grid()
plt.legend()
ax2 = ax1.twinx()
angles = np.unwrap(np.angle(w))
plt.plot(h, angles, 'g')
plt.ylabel('Angle (radians)', color='g')
plt.show()