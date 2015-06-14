#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

# Simple data to display in various forms
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

plt.close('all')

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(nrows=4, ncols=4)
# fig, axes = plt.subplots(nrows=4, ncols=4)
axarr[0].plot(x, y)
axarr[0].set_title('Sharing X axis')
axarr[1].scatter(x, y)

# # Two subplots, unpack the axes array immediately
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# ax1.plot(x, y)
# ax1.set_title('Sharing Y axis')
# ax2.scatter(x, y)

# # Three subplots sharing both x/y axes
# f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
# ax1.plot(x, y)
# ax1.set_title('Sharing both axes')
# ax2.scatter(x, y)
# ax3.scatter(x, 2 * y ** 2 - 1, color='r')
# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# f.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(18, 8))
# fig.tight_layout() # Or equivalently,  "plt.tight_layout()"

# plt.show()
# plt.figure(num=None, figsize=(18, 8), dpi=80, facecolor='w', edgecolor='k')

plt.show()