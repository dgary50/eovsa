# Change log
# 2014-06-15  DG
#   Removed old dependence on pyslalib, and added get() method 'tstamp' to
#   set or return time as a LabVIEW timestamp (seconds since 1904-01-01).  Also
#   corrected some problems, such as datime() not returning the microseconds
#   portion of the current time.
 
import datetime as dt
from time import gmtime

class datime():
    """Extend the datetime class to add mjd (Modified Julian Day) handling and
        some other conveniences.  I wanted to inherit from datetime, but that
        did not work because the datetime object, once created, is immutable,
        so I make a new date class with a datetime object dt as data, and mjd
        is calculated.
       
        Some examples: 
            d = datime()        sets d.dt to the current date and time
            d = datime(datetime(year, month[, day[, hour[, minute[, second[, us]]]]])) sets d.dt to given date/time 
            mjd = d.get()       returns the mjd (including fraction of a day)
            d.set(mjd)          sets the datetime in d.dt to the date corresponding to mjd
            d.set('2012-03-18 3:27:34.123','str')   sets the datetime in d.dt to the date in the string
            d.set('3:27:34.123','str')   sets the datetime in d.dt to today's date and the time in the string
            d.set('3:27','str')   sets the datetime in d.dt to today's date and HH:MM in the string
            st = d.get('str')   sets st to a string YYYY-MM-DD HH:MM:SS.SSSSSS
            st = d.get('mstr')  sets st to a string YYYY-Mon-DD HH:MM:SS.SSSSSS
            st = d.get('tstamp') sets st to a LabVIEW timestamp (seconds since 1904-01-01)
        
        Of course, all of the methods of datetime are available for date.dt also.
        """

    global SECONDSPERDAY

    SECONDSPERDAY = 86400.

    def __init__(self, dtin=None):
        """On init, use the passed value of a datetime object, or now() to
            populate mjd and possibly correct for time zone so that the datime
            object is in UTC.
            """
        # Check whether object was created with default "now" (done by comparing with
        # new "now", true if within 1 msec)
        if dtin is None:
            dtin = dt.datetime.now()
            
        now = dt.datetime.now()
        if (now - dtin).microseconds < 1000:
            # Was called with default "now", so compare hour with Greenwich Mean Time
            # to determine time zone
            tz = gmtime().tm_hour - dtin.hour
            # Make sure time zone is within -12 to 12 hours
            if tz < -12:
                tz = tz + 24
            elif tz > 12:
                tz = tz - 24
        else:
            tz = 0   # Do not consider time zone correction.

        # Attach datetime() object to self
        self.dt = dtin
        # Calculate mjd day using slalib routin caldj, and add time as fraction of a day
        tsec = self.dt.hour*3600. + self.dt.minute*60. + self.dt.second + self.dt.microsecond*1.e-6
        tsec
        blah = self.dt - dt.datetime(1858,11,17,0,0,0)
        mjd = blah.days + blah.seconds/86400. + blah.microseconds/86400000000.
        # Set our mjd (first call creates 
        self.__mjd = mjd
        # Now correct for time zone so that time is in UTC.
        self.set(mjd + tz/24.)

    def get(self,units='mjd'):
        """Return the date and time in various formats specified by units string
            'mjd' -- Modified Julian Day number, including fraction of a day
            'str' -- standard datetime output string, YYYY-MM-DD HH:MM:SS.SSSSSS
            'mstr'  -- string with 3-char month text,  YYYY-Mon-DD HH:MM:SS.SSSSSS
            'trad'  -- time only, as an angle in radians
            'tstamp' -- return time as LabVIEW Timestamp"""
            
        if units == 'mjd':
            """Return the Modified Julian Day number, including fraction of a day"""
            return self.__mjd
        elif units == 'str':
            """Output a standard Date/Time string (ms resolution, which can be truncated 
                by using slices if only date, or less time resolution is desired."""
            # Output string should have time to ms resolution, and the '.000' and [:23] is arranged to
            # guarantee this.
            return (self.dt.__str__() + '.000')[:23]
        elif units == 'mstr':
            """Output standard Date/TIme string except with 3-char month abbreviation"""
            return (self.dt.strftime('%Y-%b-%d %H:%M:%S.%f') + '.000')[:24]
        elif units == 'tstamp':
            """Output LabVIEW timestamp (seconds since 1904-01-01)"""
            return (self.__mjd - 16480.0)*86400.

    def set(self, value, units='mjd'):
        """Set the date and time, with input in various formats specified by units string
            'mjd' -- input is Modified Julian Day number, including fraction of a day
            'str' -- standard datetime output string, YYYY-MM-DD HH:MM:SS.SSSSSS
            'mstr'  -- string with 3-char month text,  YYYY-Mon-DD HH:MM:SS.SSSSSS
            'trad'  -- time only, as an angle in radians
            'tstamp' -- time is a LabVIEW timestamp"""
            
        if units == 'mjd':
            """Convert mjd to date/time tuple, then convert
                decimal fraction of day to HH MM SS.SSSSSS"""
            imjd = int(value)
            fmjd = (value - imjd)*86400.
            self.dt = dt.datetime(1858,11,17,0,0,0)+dt.timedelta(imjd,fmjd)
            self.__mjd = value
        elif units == 'str':
            """Use StringUtil routine secFromDMSStr() to convert time to fraction
                of a day, and convert date to days.  Value string inputs are one of
                'YYYY-MM-DD HH:MM:SS.SSS' or 'HH:MM:SS.SSS', the latter assumed to be
                for today."""
            try:
                v_date, v_time = value.split(' ')
            except:
                # Failure indicates that no date was given, so get today's date
                d = dt.datetime.today().__str__()
                v_date = d.split(' ')[0]
                v_time = value
                
            yr, mo, da = v_date.split('-')
            try:
                hh, mm, ss = v_time.split(':')
            except:
                # Failure indicates that no seconds were given, so handle
                # hour and minute, setting seconds to zero.
                hh, mm = v_time.split(':')
                ss = '0'
            try:
                ints, ms = ss.split('.')
                us = int(ms)*1000
            except:
                # Failure indicates that no fractional seconds were given,
                # so set microseconds to zero.
                ints = ss
                us = 0
            # Create a datetime object for this time, and initialize self with it.  Adding 1
            # us appears to be necessary to keep a time, e.g. 19:32:49.000000 from going to
            # 19:32:48.999999.
            self.__init__(dt.datetime(int(yr),int(mo),int(da),int(hh),int(mm),int(ints),us+1))
        elif units == 'tstamp':
            mjd = value/86400. + 16480.0  # Convert LabVIEW timestamp to mjd
            imjd = int(mjd)
            fmjd = (mjd - imjd)*86400.
            self.dt = dt.datetime(1858,11,17,0,0,0)+dt.timedelta(imjd,fmjd)
            self.__mjd = mjd
            

