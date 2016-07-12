#!/usr/bin/env python
from __future__ import division

__author__ = 'rtaborda'

import numpy as np
import matplotlib.pyplot as plt
import math
from stools import correct_baseline, scale_signal, integrate, s_filter

class seism_signal(object):
    """
    This class implements attributes and methods related to
    a single seismic record. Here, a record is understood as
    a signal or seismogram, which can be any of three types:
    an acceleration, velocity, or displacement record with
    a given sampling rate size (delta t) and a total number
    of samples
    """

    # Complete Data Group = 3-column numpy array for displ, velo, accel
    record_type = {'a': 'Acceleration', 'v': 'Velocity',
                   'd': 'Displacement', 'c': 'Complete Data Group'}

    def set_samples(self, samples):
        if not isinstance(samples, int):
            print "\nError with samples type.\n"
        self.samples = samples
    #end set_samples

    def set_dt(self, delta_t):
        if not isinstance(delta_t, float):
            print "\nError with time-step (dt).\n"
        self.dt = delta_t
    #end set_dt

    def set_data(self, data):
        # check if the data passed is a numpy array
        if not isinstance(data, np.ndarray):
            print "\nError with signal data: not a numpy array.\n"
        self.data = data
    #end set_data

    def set_type(self, stype):
        if not isinstance(stype, str) or stype not in self.record_type:
            print "\nError with signal type (must be: a, v, or d).\n"
        self.type = stype
    #end set_type

    def __init__(self, *args, **kwargs):
        """
        Initialize the record structure with all attributes empty. Can be
        initialize with unlabeled and labeld arguments Correct order for
        unlabeled arguments is: samples, dt, data, and signal type If no
        parameters, samples set to 0, dt set to 0.01, type is accelerogram and
        array is empty Acceptable
        """
        # Initialize to default values
        self.samples = 0
        self.dt = 0.01
        self.data = np.array([], float)
        self.type = 'a'

        if len(args) > 0:
            args_range = range(len(args))
            if 0 in args_range:
                self.set_samples(args[0])
            if 1 in args_range:
                self.set_dt(args[1])
            if 2 in args_range:
                self.set_data(args[2])
            if 3 in args_range:
                self.set_type(args[3])
                # all arguments were given in unlabled format
                # return

        if len(kwargs) > 0:
            if 'samples' in kwargs:
                self.set_samples(kwargs['samples'])
            if 'dt' in kwargs:
                self.set_dt(kwargs['dt'])
            if 'data' in kwargs:
                self.set_data(kwargs['data'])
            if 'type' in kwargs:
                self.set_type(kwargs['type'])

        return
    #end __init__

    def __repr__(self):
        return "Seismic record: " + \
            self.record_type[self.type] + "\n" + \
            "Samples: %i" % self.samples + "\n" + \
            "Delta t: %.3f" % self.dt + "\n" + \
            "Data:\n" + '\n'.join(map(str, self.data))
    #end __repr__

    def plot(self, flag):
        t = np.arange(0, self.samples*self.dt, self.dt)
        plt.plot(t, self.data)
        if flag == 's':
            plt.show()
        elif flag == 'p':
            # TODO: print vector file
            pass
        return
    #end plot
#end signal class

class seism_record(seism_signal):
    """
    This class extends the signal class to have addtitional
    attributes regarding time stamp and orientation
    """

    def __init__(self, *args, **kwargs):
        """
        Correct order for unlabeled arguments is: samples, dt, data, signal
        type, station name, location_lati, location_longi, depth, date, time,
        orientation
        """
        super(seism_record, self).__init__(*args, **kwargs)

        # initialize instances
        self.station_name = ""
        self.location_lati = ""
        self.location_longi = ""
        self.depth = 0.0
        self.date = ""
        self.time = ""
        self.hour = 0.0
        self.minute = 0.0
        self.seconds = 0.0
        self.fraction = 0.0
        self.tzone = ""
        self.orientation = ""

        if len(args) > 0:
            args_range = range(len(args))
            # setting location
            if 4 in args_range:
                self.set_station(args[4])
            if 5 in args_range:
                self.set_latitude(args[5])
            if 6 in args_range:
                self.set_longitude(args[6])
            if 7 in args_range:
                self.set_depth(args[7])
            if 8 in args_range:
                self.set_date(args[8])
            if 9 in args_range:
                self.set_time(args[9])
                self.set_tstamp(self.time)
            if 10 in args_range:
                self.set_orientation(args[10])
                # return
                # all arguments were given in unlabled format

        if len(kwargs) > 0:
            if 'station' in kwargs:
                self.set_station(kwargs['station'])
            if 'latitude' in kwargs:
                self.set_latitude(kwargs['latitude'])
            if 'longitude' in kwargs:
                self.set_longitude(kwargs['longitude'])
            if 'depth' in kwargs:
                self.set_depth(kwargs['depth'])
            if 'date' in kwargs:
                self.set_date(kwargs['date'])
            if 'time' in kwargs:
                self.set_time(kwargs['time'])
                self.set_tstamp(self.time)
            if 'orientation' in kwargs:
                self.set_orientation(kwargs['orientation'])

        return
    #end __init__

    def set_station(self, station_name):
        if not isinstance(station_name, str):
            print "\n**Error with station name.**\n"
        self.station_name = station_name
    #end set_station

    def set_latitude(self, location_lati):
        # checking latitude format being float+S/N
        if not isinstance(location_lati, str):
            print "\n**Error with location latitude (Invalid instance type).**\n"
        elif location_lati[-1] not in ["N", "S"]:
            print "\nError with location latitude (Invalid format).\n"
        else:
            try:
                float(location_lati[0:-2])
            except ValueError:
                print "\nError with location latitude (Invalid format).\n"
        self.location_lati = location_lati
    #end set_latitude

    def set_longitude(self, location_longi):
        # checking longitude format being float+E/W
        if not isinstance(location_longi, str):
            print "\n**Error with location longitude (Invalid instance type).**\n"
        elif location_longi[-1] not in ["E", "W"]:
            print "\n**Error with location longitude (Invalid format).**\n"
        else:
            try:
                float(location_longi[0:-2])
            except ValueError:
                print "\n**Error with location longitude (Invalid format).**\n"

        self.location_longi = location_longi
    #end set_longitude

    def set_depth(self, depth):
        # checking depth
        if not isinstance(depth, float):
            print "\n**Error with depth.**\n"
        self.depth = depth
    #end set_depth

    def set_orientation(self, orientation):
        # if the orientation is string, it should be either Up or Down
        # if the orientation is int, it should between 0 and 360
        if isinstance(orientation, str) and orientation.upper() in ["UP", "DOWN"]:
            self.orientation = orientation
        elif isinstance(orientation, int) and orientation <= 360 and orientation >= 0:
            self.orientation = orientation
        else:
            print "[ERROR]: Invalid orientation."
    #end set_orientation

    def set_date(self, date):
        # check the format of date string being #/#/#
        if not isinstance(date, str):
            print "[ERROR]: invalid date."
        else:
            for x in date.split('/'):
                if x.isdigit() == False:
                    print "[ERROR]: invalid date."
                    break
        self.date = date

    def set_time(self, time):
        if not isinstance(time, str):
            print "[ERROR]: invalid date."
        self.time = time

    # the function is to split time string into hour, minute, seconds, fraction, and tzone
    def set_tstamp(self, time):
        hour = self.time.split(':')[0]
        minute = self.time.split(':')[1]
        seconds = time.split(':')[2].split()[0].split('.')[0]
        fraction = time.split(':')[2].split()[0].split('.')[1]
        tzone = self.time.split()[-1]

        try:
            hour = float(hour)
            minute = float(minute)
            seconds = float(seconds)
            fraction = float(fraction)
        except ValueError:
            print "[ERROR]: invalid start time."

        self.hour = hour
        self.minute = minute
        self.seconds = seconds
        self.fraction = fraction
        self.tzone = tzone

    # to test with record object
    def print_attr(self):
        print "==================================================================="
        print "samples: " + str(self.samples)
        print "dt: " + str(self.dt)
        print "data type: " + self.type
        print self.data
        print "station name: " + self.station_name
        print "station latitude: " + self.location_lati
        print "station longitude: " + self.location_longi
        print "depth: ??"
        print "date: " + self.date
        print "time: " + self.time
        print "hour: " + str(self.hour)
        print "minute: " + str(self.minute)
        print "seconds: " + str(self.seconds)
        print "fraction: " + str(self.fraction)
        print "tzone: " + self.tzone
        print "orientation: " + str(self.orientation)
    #end print_attr

    def process_ori(self):
        # process data by orientation
        if isinstance(self.orientation, int):
            if self.orientation in [0, 360, 90, -270]:
                pass
            elif self.orientation in [180, -180, -90, 270]:
                # if negative/down orientation, multiply by -1
                self.data = self.data*-1
                self.orientation += 180

                if self.orientation > 360:
                    self.orientation -= 360

            else:
                # encounter special orientation; pass to rotate
                return True

        elif isinstance(self.orientation, str):
            if self.orientation == 'Up':
                pass
            elif self.orientation == 'Down':
                self.data = self.data*-1
                self.orientation = 'Up'
            else:
                # encounter invalid orientation
                return False
    # end of process_ori

    # def convert(self):
    #     # make average on first 10% of samples; minus average and multiply by 981
    #     self.data = 981*(self.data - np.average(self.data[0:int(self.samples*0.1)]))
    #     pass

    # def process_smc_v1(self):
    #     """
    #     The function process record by converting orientation and multiplying data
    #     """
    #     # If the orientation was not set properly, it would be empty string by default
    #     # if not self.orientation:
    #     #     pass
    #         # print "[ERROR]: invalid orientation"
    #         # return False
    #         # orientation = " "
    #     if self.orientation in [0, 360]:
    #         orientation = 'N'
    #     elif self.orientation in [180, -180]:
    #         orientation = 'N'
    #         self.data = self.data*-1
    #     elif self.orientation in [90, -270]:
    #         orientation = 'E'
    #     elif self.orientation in [-90, 270]:
    #         orientation = 'E'
    #         self.data = self.data*-1
    #     elif self.orientation == "Up":
    #         orientation = 'Z'
    #     elif self.orientation == "Down":
    #         orientation = 'Z'
    #         self.data = self.data*-1
    #     elif isinstance(self.orientation, int):
    #         # orientation needed to be roated
    #          return True
    #     else:
    #         # invalid orientation
    #         return False

    # make average on first 10% of samples; minus average and multiply by 981
    # self.data = 981*(self.data - np.average(self.data))
    # self.data = 981*(self.data - np.average(self.data[0:int(self.samples*0.1)]))
    # return
    #end process_smc_V1

# end record class

# ================ Classes for processed data (V2) =================
class seism_psignal(seism_signal):
    """
    This class extends the seism_signal class to have addtitional
    attributes for acceleration, velocity, and displacement data
    """

    def __init__(self, *args, **kwargs):
        """
        Correct order for unlabeled arguments is: samples, dt, data,
        signal type (eliminate), acceleration, velocity, displacement
        """
        super(seism_psignal, self).__init__(*args, **kwargs)

        # initialize accel, velo, displ
        self.accel = np.array([], float)
        self.velo = np.array([], float)
        self.displ = np.array([], float)


        if len(args) > 0:
            args_range = range(len(args))
            if 4 in args_range:
                self.set_accel(args[4])
            if 5 in args_range:
                self.set_velo(args[5])
            if 6 in args_range:
                self.set_displ(args[6])
                return

        if len(kwargs) > 0:
            if 'accel' in kwargs:
                self.set_accel(kwargs['accel'])
            if 'velo' in kwargs:
                self.set_velo(kwargs['velo'])
            if 'displ' in kwargs:
                self.set_displ(kwargs['displ'])
        return
    #end __init__

    def set_accel(self, accel):
        # check if the data passed is a numpy array
        if not isinstance(accel, np.ndarray):
            print "[ERROR]: signal acceleration data - not a numpy array."
        self.accel = accel
    #end set_accel

    def set_velo(self, velo):
        if not isinstance(velo, np.ndarray):
            print "[ERROR]: signal velocity data - not a numpy array."
        self.velo = velo
    #end set_velo

    def set_displ(self, displ):
        if not isinstance(displ, np.ndarray):
            print "[ERROR]: signal displacement data - not a numpy array."
        self.displ = displ
    #end set_displ

    def print_attr(self):
        print "===========================psignal========================================"
        print "samples: " + str(self.samples)
        print "dt: " + str(self.dt)
        print "data type: " + self.type
        print self.data
        print self.accel
        print self.velo
        print self.displ

#end seism_psignal class

class seism_precord(seism_record):
    """
    This class extends the seism_record class to have addtitional
    attributes for acceleration, velocity, and displacement data
    """

    def __init__(self, *args, **kwargs):
        """
        Correct order for unlabeled arguments is: samples, dt, data, signal type
        (eliminate), station, location_lati, location_longi, depth, date, time,
        orientation, acceleration, velocity, displacement
        """
        super(seism_precord, self).__init__(*args, **kwargs)

        # initialize accel, velo, displ
        self.accel = np.array([], float)
        self.velo = np.array([], float)
        self.displ = np.array([], float)

        if len(args) > 0:
            args_range = range(len(args))
            if 11 in args_range:
                self.set_accel(args[11])
            if 12 in args_range:
                self.set_velo(args[12])
            if 13 in args_range:
                self.set_displ(args[13])
                return

        if len(kwargs) > 0:
            if 'accel' in kwargs:
                self.set_accel(kwargs['accel'])
            if 'velo' in kwargs:
                self.set_velo(kwargs['velo'])
            if 'displ' in kwargs:
                self.set_displ(kwargs['displ'])
        return
    #end __init__

    def set_accel(self, accel):
        # check if the data passed is a numpy array
        if not isinstance(accel, np.ndarray):
            print "\n[ERROR]: signal acceleration data - not an numpy array.\n"
        self.accel = accel
    #end set_accel

    def set_velo(self, velo):
        if not isinstance(velo, np.ndarray):
            print "\n[ERROR]: signal velocity data - not an numpy array.\n"
        self.velo = velo
    #end set_velo

    def set_displ(self, displ):
        if not isinstance(displ, np.ndarray):
            print "\n[ERROR]: signal displacement data - not an numpy array.\n"
        self.displ = displ
    #end set_displ

    def print_attr(self):
        print "=============================precord======================================"
        print "samples: " + str(self.samples)
        print "dt: " + str(self.dt)
        print "data type: " + self.type
        print "station name: " + self.station_name
        print "station latitude: " + self.location_lati
        print "station longitude: " + self.location_longi
        print "depth: ??"
        print "date: " + self.date
        print "time: " + self.time
        print "hour: " + str(self.hour)
        print "minute: " + str(self.minute)
        print "seconds: " + str(self.seconds)
        print "fraction: " + str(self.fraction)
        print "tzone: " + self.tzone
        print "orientation: " + str(self.orientation)
        print self.accel
        print self.velo
        print self.displ
#end seism_precord class

# ==========================Class for Station==============================
class seism_station(object):
    """
    The station object contains a list of signals/psignals;
    with network number and station ID
    """
    def __init__(self, *args, **kwargs):
        # Initialize to default values
        self.list = []
        self.network = ''
        self.id = ''
        self.type = ''
        self.name = ''
        self.latitude = ''
        self.longitude = ''


        if len(args) > 0:
            args_range = range(len(args))
            if 0 in args_range:
                if self.set_list(args[0]) == False:
                    return
            if 1 in args_range:
                self.set_network(args[1])
            if 2 in args_range:
                self.set_id(args[2])
            if 3 in args_range:
                self.set_type(args[3])
                # all arguments were given in unlabled format
                # return

        if len(kwargs) > 0:
            if 'list' in kwargs:
                if self.set_list(kwargs['list']) == False:
                    return
            if 'network' in kwargs:
                self.set_network(kwargs['network'])
            if 'station' in kwargs:
                self.set_id(kwargs['station'])
            if 'type' in kwargs:
                self.set_type(kwargs['type'])

        self.set_name()
        self.set_location()

        return
    #end __init__

    def set_list(self, record_list):
        if len(record_list) != 3:
            print "[ERROR]: The program handles stations with three channels ONLY. "
            return False

        if not all(isinstance(record, seism_signal) for record in record_list):
            print "[ERROR]: list in station does not contain seism signals."
            return False

        if all(isinstance(record, seism_record) for record in record_list):
            if not (record_list[0].orientation
                    != record_list[1].orientation
                    != record_list[2].orientation):
                print "[ERROR]: conflict orientations."
                return False

        self.list = record_list
    # end of set_list

    def set_network(self, network):
        self.network = network

    def set_id(self, station_id):
        self.id = station_id

    def set_type(self, filetype):
        filetype = filetype.upper()
        if filetype in ['V1', 'V2']:
            self.type = filetype

    def set_name(self):
        if self.list:
            if (self.list[0].station_name \
                == self.list[1].station_name \
                == self.list[2].station_name):
                self.name = self.list[0].station_name
            else:
                print "[ERROR]: channels have different station names."
                return
        return
    # end of set_name

    def set_location(self):
        if self.list:
            if (self.list[0].location_lati \
                == self.list[1].location_lati \
                == self.list[2].location_lati):
                self.latitude = self.list[0].location_lati

            if (self.list[0].location_longi \
                == self.list[1].location_longi \
                == self.list[2].location_longi):
                self.longitude = self.list[0].location_longi

            else:
                print("[ERROR]: channels have different "
                      "location (latitude/longitude).")
                return
        return
    # end of set_location

    def rotate(self, record_list, flag):
        """
        The function is to transfrom data for
        channels with special orientations.
        """
        tmp = []
        for record in record_list:
            if record.orientation != 'Up':
                tmp.append(record)
                record_list.remove(record)
        if len(tmp) != 2:
            # Couldn't find two channels to rotate
            return False
        x = tmp[0].orientation
        y = tmp[1].orientation

        if abs(x-y) != 90:
            # We need two orthogonal channels
            return False

        if x > y:
            list(reversed(tmp))

        if flag == 'v1':
            # rotate
            matrix = np.array([(math.cos(math.radians(x)),
                                -math.sin(math.radians(x))),
                               (math.sin(math.radians(x)),
                                math.cos(math.radians(x)))])
            data = matrix.dot([tmp[0].data, tmp[1].data])

            # transform the first record with North orientation
            tmp[0].data = data[0]
            tmp[0].orientation = 0

            # transform the second record with East orientation
            tmp[1].data = data[1]
            tmp[1].orientation = 90

            record_list += tmp
            return record_list

        if flag == 'v2':
            matrix = np.array([(math.cos(math.radians(x)),
                                -math.sin(math.radians(x))),
                               (math.sin(math.radians(x)),
                                math.cos(math.radians(x)))])
            [tmp[0].accel, tmp[1].accel] = matrix.dot([tmp[0].accel,
                                                       tmp[1].accel])
            [tmp[0].velo, tmp[1].velo] = matrix.dot([tmp[0].velo,
                                                     tmp[1].velo])
            [tmp[0].disp, tmp[1].disp] = matrix.dot([tmp[0].displ,
                                                     tmp[1].displ])

            tmp[0].orientation = 0
            tmp[1].orientation = 90

            record_list += tmp
            return record_list

    # end of rotate

    def process_v1(self):
        """
        The function is take a list of records get from V1 files;
        then use their data to all three types of data (acc, vel, dis),
        then return a list of precords.
        """
        rotate_flag = False
        precord_list = []

        for record in self.list:
            # process data of record
            # if not record.orientation:
            #     print record.orientation
            #     print "[ERROR]: invalid orientation."
            #     return False


            # reverse the data by orientation
            if record.process_ori() == True:
                # if encounter special orientations.
                rotate_flag = True
                break
            elif record.process_ori() == False:
                # if encouter invalid orientations.
                return False

            # record.convert() # convert the unit of data
            correct_baseline(record)
            scale_signal(record, 981)

            if record.type == 'a':
                # get velocity and displacement
                record.data = s_filter(record.data, record.dt,
                                       type='highpass', family='butter',
                                       fmin=0.05, N=5)

                velocity = integrate(record.data, record.dt)
                velocity = s_filter(velocity, record.dt,
                                    type='highpass', family='butter',
                                    fmin=0.05, N=5)

                displacement = integrate(velocity, record.dt)
                displacement = s_filter(displacement, record.dt,
                                        type='highpass', family='butter',
                                        fmin=0.05, N=5)

                data = np.c_[displacement, velocity, record.data]

                precord = seism_precord(record.samples, record.dt, data,
                                        'c', record.station_name,
                                        accel=record.data,
                                        displ=displacement, velo=velocity,
                                        orientation=record.orientation,
                                        date=record.date, time=record.time,
                                        depth=record.depth,
                                        latitude=record.location_lati,
                                        longitude=record.location_longi)
                precord_list.append(precord)

        # rotation
        if rotate_flag:
            record_list = self.rotate(self.list, 'v1')
            if not record_list:
                return False
            else:
                self.list = record_list
            # recursively calling the function to continue processing
            return self.process_v1()

        self.list = precord_list
        return True
    # end of process_v1

    def process_v2(self):
        """
        Checks records from V2 files; rotate if necessary so the
        horizontal channels end up N/S and E/W
        """
        rotate_flag = False
        for record in self.list:
            if (isinstance(record.orientation, int)) \
                and (record.orientation not in [0, 360, 180, -180,
                                                90, 270, -90, -270]):
                rotate_flag = True
                break

        if rotate_flag:
            record_list = self.rotate(self.list, 'v2')
            if not record_list:
                return False
            else:
                self.list = record_list
                return True
        else:
            # rotation not needed
            return True
    # end of process_v2

#end station class
