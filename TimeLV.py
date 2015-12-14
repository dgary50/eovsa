#  Extend the astropy Time object to work with LabVIEW times.
#    usage: from Time_LV import Time_LV
#           t = Time_LV([3513711600.234,3513715140.359],0.001,format='lv')
#           print t.lv
#           [  3.51371160e+09   3.51371514e+09]
#           print t.iso
#           ['2015-05-05 23:00:00.235' '2015-05-05 23:59:00.360']
#  2015-05-24  DG
#    First written.
#  2015-05-26  DG
#    Change class name from Time_LV to just Time, in analogy to the
#    astropy.time Time.
#

from astropy.time import Time as astroTime
import numpy as np

from astropy.time.core import TimeFromEpoch
class TimeLV(TimeFromEpoch):
    """
    LabVIEW timestamp input:
    number of seconds from 1904-01-01 00:00:00 UTC

    """
    # This corresponds to the zero reference time for LabVIEW.
    # Note that TAI and UTC are equivalent at the reference time.
    name = 'lv'
    unit = 1.0/86400.0     # Converts to seconds
    epoch_val = 2416480.5  # Time('1904-01-01 00:00:00', scale='tai').jd - 1
    epoch_val2 = None
    epoch_scale = 'utc'
    epoch_format = 'jd'
    
class Time(astroTime):
    ''' Extends astropy Time object to handle LabVIEW timestamps
        (format 'lv').
    '''
    def __init__(self, val, val2=None, format=None, scale=None,
             precision=None, in_subfmt=None, out_subfmt=None,
             location=None, copy=False):
        # Extend formats list to include TimeLV (LabVIEW format)
        self.FORMATS[u'lv'] = TimeLV
        astroTime.__init__(self, val, val2, format=format, scale=scale,
             precision=precision, in_subfmt=in_subfmt, out_subfmt=out_subfmt,
             location=location, copy=copy)
