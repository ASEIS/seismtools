#!/usr/bin/env python
import os 
print os.stat("SampleFiles/foo").st_size == 0