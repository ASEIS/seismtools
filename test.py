#!/usr/bin/env python
column1 = ["soft","pregnant","tall"]
column2 = ["skin","woman", "man"]

for c1, c2 in zip(column1, column2):
    print "%-9s %s" % (c1, c2)