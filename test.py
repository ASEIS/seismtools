#!/usr/bin/env python
import numpy as np
import math
y = 104 
x = 14
x = (x+y-90)/2 
h1 = 20 
h2 = 90

# x = 104 
# y = 14
# h1 = 20 
# h2 = 90


N = h1*math.cos(math.radians(x)) - h2*math.sin(math.radians(x))
E = h1*math.sin(math.radians(x)) + h2*math.cos(math.radians(x))
print N 
print E


j = np.array([(math.cos(math.radians(x)), -math.sin(math.radians(x))), (math.sin(math.radians(x)), math.cos(math.radians(x)))])
# w = np.array([h1, h2])
print(j.dot([h1, h2]))
print(j.dot([h2, h1]))

