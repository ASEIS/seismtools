
from seism import *

def load_smc_v1(filename):
    
    # loads station into a string
    fp = open(filename, 'r')
    channels = fp.read()
    fp.close()

    # splits the string by channels
    channels = channels.split('/&')
    del(channels[len(channels)-1])

    # splits the channels
    for i in range(len(channels)):
        channels[i] = channels[i].split('\r\n')

    # clean the first row in all but the first channel
    # this row corresponds to the channel delimiter
    for i in range(1,len(channels)):
        del channels[i][0]

    for i in range(len(channels)):

        # check this is the uncorrected acceleration data
        ctype = channels[i][0][0:24].lower()
        if ctype != "uncorrected accelerogram":
            return channels, 0

        tmp = channels[i][27].split()
        samples = int(tmp[0])

        tmp = channels[i][28:]
        signal = str()
        for s in tmp:
            signal += s
        signal = signal.split()

        record = seism_record(samples=samples, dt=0.005, data=signal, type='a')
        print record
        pass

    return channels, 1

