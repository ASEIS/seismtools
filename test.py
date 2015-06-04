#!/usr/bin/env python
try:
    fp = open("nother")
except IOError, e:
    # print e.errno
    print e

