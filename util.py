# 2014-11-17  Operator  <sched@helios>
#   Made changes to datime to allow setting/getting time as LabVIEW 
#   timestamp.
# 2015-01-09  DG
#   Added 1-s time-out for reading IERS bulletin
# 2015-05-23  DG
#   Changed datime() to return mjd on set (used to return nothing)
# 2015-05-27  DG
#   Added a new Time class, that extends the astropy Time object to 
#   work with LabVIEW times.
#    usage: from util import Time
#           t = Time([3513711600.234,3513715140.359],0.001,format='lv')
#           print t.lv
#           [  3.51371160e+09   3.51371514e+09]
#           print t.iso
#           ['2015-05-05 23:00:00.235' '2015-05-05 23:59:00.360']
#  2015-Jun-28  DG
#    Added common_val_idx() from solpnt, so that I do not have to include
#    the entire solpnt just for this!
#  2016-Oct-19  DG
#    Added nearest_val_idx(), which looks in second array for nearest values
#    to those in first array.
#  2017-Jan-15  DG
#    Added lobe() function to put phase into +/- pi range
#  2017-Aug-08  DG
#    Added precision keyword to common_val_idx()
#    Also added code for calculating baseline order
#  2019-Jan-03  DG
#    Added lin_phase_fit() routine, which is a common need to find the
#    phase slope of a phase spectrum.
#  2019-Feb-18  DG
#    Added sat_phase, which determines corrections for XY and YX phase based
#    on a packet dump of an R/L polarized communications satellite (e.g. Ciel)
#  2019-Jul-15  DG
#    Added plot_sched_log() routine to plot the times and durations of non-standard
#    durations in the /tmp/schedule.log file.
#  2019-Jul-18  DG
#    Added fname2mjd(), versions of which have been in several other modules.
#    This version allows conversion of arrays or lists of filenames.
#  2019-Aug-24  DG
#    Fixed a bug in lin_phase_fit() when data were all nan.
#  2020-Jan-30  DG
#    Fixed another bug in lin_phase_fit() that caused shape error.
#  2020-May-02  DG
#    Added fix_time_drift() routine based on the one with the same name in
#    calwidget.py (although that one assumes a different data structure). This
#    verstion works with a standard read_idb() output file.
#  2020-May-10  DG
#    Added routine get_idbdir(), which returns the root path of IDB data for a
#    given date.
#  2020-May-10  SY
#    Added envirmental variable EOVSADBJSON to get_idbdir().
#    Use the json file defined by EOVSADBJSON to control the date-dependant data location.
#  2020-May-11  DG
#    Fixed failure of get_idbdir() when EOVSADBJSON is not defined.
#  2020-May-23  DG
#    get_idbdir() was returning a uencoded string, which killed aipy. Now
#    converted to an ASCII string.
#  2020-Aug-02  DG
#    Apparently the IERS bulletin URL changed some time ago and finally today
#    the old file no longer had information for today's date (the stale file
#    from about October of last year was being used since then!).  I found the
#    new URL and updated the UT1-UTC() function.  No other change was needed.
#  2021-Jan-01  DG
#    Added freq2bdname(), to call function of the same name from either chan_util_bc
#    or chan_util_52, depending on the date.
#  2022-Jun-23  DG
#    Added suppression of spurious ERFA warnings involving absence of leap
#    seconds for some (mainly future) dates.
# *

import StringUtil as su
from numpy import pi, sqrt, array, mat, matrix, dot, where, ndarray
import datetime as dt
from time import gmtime
#
# Suppress spurious warnings from ERFA
import warnings
from astropy._erfa.core import ErfaWarning
warnings.simplefilter('ignore', category=ErfaWarning)

class Angle:
    """General angle class, converts input to radians, but handles input
       and output in various units."""

    __pi = pi
    __twopi = 2 * pi
    __degperrad = 180. / pi
    __arcsecperdeg = 3600.
    __masperdeg = 3600000.
    __arcsecperrad = __arcsecperdeg * __degperrad
    __arcminperrad = 60 * __degperrad
    __masperrad = __masperdeg * __degperrad

    def __init__(self, value=0.0, units='radians'):
        """Constructor"""
        self.modulo = False
        self.set(value, units)

    def set(self, value=0.0, units='radians'):
        """Sets radians value, converting according to units string
            'radians' -- no conversion
            'degrees' -- decimal degrees to radians
            'dmsstr'  -- sexagesimal string degress as -DD:MM:SS.SSS to radians
            'hmsstr'  -- sexagesimal string hours as   -HH:MM:SS.SSS to radians
            'masec'   -- milliarcseconds to radians
            'arcsec'  -- arcseconds to radians
            'dms' and 'hms' are aliases for 'dmsstr' and 'hmsstr' """

        if units == 'radians':
            self.radians = value
        elif units == 'degrees':
            self.radians = value / self.__degperrad
        elif units == 'dmsstr' or units == 'dms':
            # String in the form -DD:MM:SS.SSS
            self.radians = su.secFromDMSStr(value) / self.__arcsecperrad
        elif units == 'hmsstr' or units == 'hms':
            # String in the form -HH:MM:SS.SSS
            self.radians = su.secFromDMSStr(value) * 15. / self.__arcsecperrad
        elif units == 'masec':
            self.radians = value / self.__masperrad
        elif units == 'arcsec':
            self.radians = value / self.__arcsecperrad
        else:
            print 'ERROR: Unknown units string, radians used'
            self.radians = value

        # If modulo bool is set, make sure the value of radians
        # is between 0, 2pi
        if self.modulo:
            radians = self.radians
            twopi = self.__twopi
            while radians >= twopi:
                radians -= twopi
            while radians < 0.0:
                radians += twopi
            self.radians = radians

    def get(self, units='radians'):
        """Gets radians value and converts according to units string
            'radians' -- no conversion
            'degrees' -- radians to decimal degrees
            'dmsstr'  -- radians to sexagesimal string degress as -DD:MM:SS.SSS
            'hmsstr'  -- radians to sexagesimal string hours as   -HH:MM:SS.SSS
            'masec'   -- radians to milliarcseconds
            'arcsec'  -- radians to arcseconds
            'dms' and 'hms' are aliases for 'dmsstr' and 'hmsstr' """

        if units == 'radians':
            return self.radians
        elif units == 'degrees':
            return self.radians * self.__degperrad
        elif units == 'dmsstr' or units == 'dms':
            # String in the form -DD:MM:SS.SSS
            strg = su.dmsStrFromDeg(self.radians * self.__degperrad, precision=3)
            if len(strg) == 11:
                # Need to add missing zero and space
                strg = ' 0' + strg
            elif len(strg) == 12:
                # Need to add missing zero or space
                if strg[0] == '-':
                    strg = strg[0] + '0' + strg[1:]
                else:
                    strg = ' ' + strg
            return strg
        elif units == 'hmsstr' or units == 'hms':
            # String in the form -HH:MM:SS.SSS
            strg = su.dmsStrFromDeg(self.radians * self.__degperrad / 15., precision=3)
            if len(strg) == 11:
                # Need to add missing zero
                strg = '0' + strg
            return strg
        elif units == 'masec':
            return self.radians * self.__masperrad
        elif units == 'arcsec':
            return self.radians * self.__arcsecperrad
        else:
            print 'ERROR: Unknown units string, radians used'
            return self.radians

    def __add__(self, other):
        """Add one Angle object to another, or add a scalar to
            an Angle object.
            """
        if isinstance(other, int) or isinstance(other, float):
            # scalar is interpreted as radians
            radians = other
        elif isinstance(other, Angle):
            radians = other.radians
        else:
            print 'Angle.__add__(): Cannot add', other, 'to an Angle()'
            return 0

        # Create the right kind of object to return (same as self)
        if isinstance(self, RA_Angle):
            return RA_Angle(self.radians + radians)
        elif isinstance(self, Dec_Angle):
            return Dec_Angle(self.radians + radians)
        else:
            return Angle(self.radians + radians)

    def __sub__(self, other):
        """Subtract one Angle object from another, or subtract a scalar from
            an Angle object.
            """
        if isinstance(other, int) or isinstance(other, float):
            # scalar is interpreted as radians
            radians = other
        elif isinstance(other, Angle):
            radians = other.radians
        else:
            print 'Angle.__sub__(): Cannot subtract', other, 'from an Angle()'
            return 0

        # Create the right kind of object to return (same as self)
        if isinstance(self, RA_Angle):
            return RA_Angle(self.radians - radians)
        elif isinstance(self, Dec_Angle):
            return Dec_Angle(self.radians - radians)
        else:
            return Angle(self.radians - radians)


class RA_Angle(Angle):
    """Inherits Angle class, with modulo = True"""

    def __init__(self, value=0.0, units='radians'):
        """Constructor"""
        self.modulo = True
        self.set(value, units)


#    def set(self, value, units='radians'):
#        Angle.set(self, value, units)


class Dec_Angle(Angle):
    """Inherits Angle class, with different modulo action"""

    def __init__(self, value=0.0, units='radians'):
        """Constructor"""
        self.modulo = True
        self.set(value, units)

    def set(self, value=0.0, units='radians'):
        # Set modulo false for a moment
        save_modulo = self.modulo
        self.modulo = False
        Angle.set(self, value, units)
        # Now set back to original (normally True unless overridden)
        self.modulo = save_modulo

        # If modulo bool is set, make sure the value of radians
        # is between -pi, pi.  This doesn't guarantee it is a
        # valid Declination, but at least gives it a fighting chance.
        if self.modulo:
            radians = self.radians
            twopi = 2. * pi
            while radians >= pi:
                radians -= twopi
            while radians < -pi:
                radians += twopi
            self.radians = radians
            # Set the decvalid flag to False if not between -pi/2 and pi/2.
            if radians >= pi / 2 or radians < -pi / 2:
                self.decvalid = False
            else:
                self.decvalid = True


class Length():
    """Class for defining lengths and their units.  The native unit for
        the object is meters.
        """

    global cmperm, mperkm

    cmperm = 100.
    mperkm = 1000.

    def __init__(self, value=0.0, units='m'):
        """Constructor allows setting a value on creation.  Default is
            a zero-length string with meter units.
            """
        self.set(value, units)

    def set(self, value=0.0, units='m'):
        """Sets the value and units of the Length object.  Default is a
            zero-length string with meter units.

            Choice of units is:
                'cm' = centimeters
                'm'  = meters
                'km' = kilometers
            """
        if not (isinstance(value, int) or isinstance(value, float)):
            print 'Length.set(): Invalid length value', value, 'Setting to zero'
            self.meters = 0.0

        if units == 'cm':
            self.meters = value / cmperm
        elif units == 'm':
            self.meters = float(value)
        elif units == 'km':
            self.meters = value * mperkm
        else:
            print 'Length.set(): Invalid length units', units, 'Meters assumed.'

        self.units = units

    def get(self, units='m'):
        """Gets the value of the Length object, in the units specified. Defaults
            to returning the value in meters, e.g.
                >>> a = Length(0.175,'km')
                >>> a.get()
                175.0

            Choice of units is:
                'cm' = return float in centimeter units
                'm'  = return float in meter units
                'km' = return float in kilometer units
                'str' = return string value in object's own units
            """
        if units == 'cm':
            value = self.meters * cmperm
        elif units == 'm':
            value = self.meters
        elif units == 'km':
            value = self.meters / mperkm

        if units == 'str':
            value = self.get(units=self.units)
            return str(value) + ' ' + self.units
        else:
            return value

    def __add__(self, other):
        """Adds two Length() objects or a Length() object and a constant in
            meter units. In either case, result has units of self.
                >>> a = Length(0.175,'km')   # Define 175 m length
                >>> b = Length(25)           # Define 25 m length
                >>> (a+b).get()
                200.0
                >>> (a + 27).get()
                202.0
            """
        if isinstance(other, Length):
            new = Length(self.meters + other.meters)
        elif isinstance(other, int) or isinstance(other, float):
            new = Length(self.meters + other)
        else:
            print 'Length.__add__(): Cannot add', other, \
                'to a Length object. Returning copy of self.'
            new = Length(self.meters)
        new.units = self.units

        return new

    def __sub__(self, other):
        """Subtracts two Length() objects 
                >>> a = Length(0.175,'km')   # Define 175 m length
                >>> b = Length(25)           # Define 25 m length
                >>> (a-b).get()              # Subtract two lengths, return meters
                150.0

            or a Length() object and a constant in meter units. 
            In either case, result has units of self.
                >>> (a - 27).get()           # Subtract constant from a length
                148.0
            """
        if isinstance(other, Length):
            new = Length(self.meters - other.meters)
        elif isinstance(other, int) or isinstance(other, float):
            new = Length(self.meters - other)
        else:
            print 'Length.__sub__(): Cannot subtract', other, \
                'from a Length object. Returning copy of self.'
            new = Length(self.meters)
        new.units = self.units

        return new

    def __div__(self, other):
        """Divides two Length() objects and returns a constant without units,
                >>> a = Length(0.175,'km')
                >>> b = Length(25)
                >>> a/b
                7.0

            or a Length() object and a constant and returns a Length() object
            in units of self.
                >>> a = Length(0.175,'km')
                >>> (a/3.1).get()
                56.45161290322581

            """
        if isinstance(other, Length):
            new = self.meters / other.meters
        elif isinstance(other, int) or isinstance(other, float):
            new = Length(self.meters / other)
            new.units = self.units
        else:
            print 'Length.__div__(): Cannot divide a Length object by', other, \
                'Error return "False"'
            new = False

        return new


class Vector():
    """A general class to hold and manipulate 3-element vectors of
        Length() objects
        """

    def __init__(self, lengths=None, units='m'):
        """Constructor, input defaults to a three-element vector of three
            zero-length Length objects with meter units
            """
        if lengths is None:
            # Handle case of initialization with no arguments
            # Believe it or not, using lengths = [0.0, 0.0, 0.0] in the argument
            # list caused the same Length() objects to be returned for every
            # Vector() call!  I have since learned that this is a "standard
            # problem." Functions in argument lists are evaluated at compile
            # time, not anew every time the function is invoked.
            lengths = [0.0, 0.0, 0.0]
        self.set(lengths, units)

    def set(self, lengths=None, units='m'):
        """Set the value of the Vector using any combination of Length()
            objects or constants, the latter interpreted according to
            the specified units.  The default is three zero-length Length()
            objects in meter units.
            """
        if lengths is None:
            lengths = [0.0, 0.0, 0.0]
        if len(lengths) != 3:
            print 'Vector.set(): Must supply three lengths. Object not set.'
            lengths = [Length(0), Length(0), Length(0)]

        # If a numpy array was passed in, change to a list
        try:
            lengths = lengths.tolist()
        except:
            pass

        # Loop over the three lengths to check them
        for i in range(3):
            if not isinstance(lengths[i], Length):
                # One of the supplied objects is not a Length() object, so
                # see if it is a scalar constant.
                if isinstance(lengths[i], int) or isinstance(lengths[i], float):
                    # It is, so create a Length() object using it and the
                    # supplied units.
                    lengths[i] = Length(lengths[i], units=units)
                else:
                    # It is neither a Length object nor a scalar constant, so
                    # it is an error.
                    print 'Vector.set(): Item', i, 'is not a Length() or a constant.\n\
                        Replacing with a zero-length Length object'
                    lengths[i] = Length()

        self.lengths = lengths
        self.units = units

    def get(self, units='m'):
        """Get the lengths in the Vector in the units specified.  Defaults
            to returning values in meters.
            """
        out = []
        # Loop over the three lengths to return their values as a list
        for i in range(3):
            # if units == 'str':
            #    value = self.lengths[i].get(units)
            #    out.append(str(value) + ' ' + self.units)
            # else:
            value = self.lengths[i].get(units)
            out.append(value)

        return out

    def __add__(self, other):
        """Add two vector objects"""
        if isinstance(other, Vector):
            out = Vector()
            for i in range(3):
                out.lengths[i] = self.lengths[i] + other.lengths[i]
            return out
        else:
            print 'Vector.__add__(): Object', other, 'must be a Vector() object.', \
                'Returning a zero vector.'
            return Vector([0.0, 0.0, 0.0])

    def __sub__(self, other):
        """Subtract two vector objects"""
        if isinstance(other, Vector):
            out = Vector()
            for i in range(3):
                out.lengths[i] = self.lengths[i] - other.lengths[i]
            return out
        else:
            print 'Vector.__sub__(): Object', other, 'must be a Vector() object.', \
                'Returning a zero vector.'
            return Vector([0.0, 0.0, 0.0])

    def magnitude(self, units='m'):
        """Returns the magnitude of the Vector in the units specified.
            """
        if units == 'str':
            # Return the lengths in the units of the Vector object
            vals = self.get(self.units)
        else:
            # Return the lengths in the specified units
            vals = self.get(units)

        # Calculate the magnitude
        mag = 0.0
        for i in range(3):
            mag += vals[i] * vals[i]
        magnitude = sqrt(mag)
        if units == 'str':
            # Make a string out of it
            return str(magnitude) + ' ' + self.units

        return magnitude

    def rotate(self, rot=None):
        """Rotate a Vector object by applying the rotation matrix rot."""

        if rot is None:
            rot = mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if not isinstance(rot, matrix):
            print 'Vector.rotate(): Object', rot, 'is not a matrix.  No rotation applied.'
            rot = mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        # Apply rotation (output is in meter)
        vnew = dot(rot, self.get())
        # Create a new Vector object for output.  Have to convert vnew numpy
        # array to a "list of lists", then take the first element
        new = Vector(vnew.tolist()[0])
        new.units = self.units  # Set output units same as input
        return new

    def __mul__(self, other):
        """Multiplies a Vector() object by a scalar.
            """
        if isinstance(other, int) or isinstance(other, float):
            new = Vector((array(self.get()) * other).tolist())
            new.units = self.units
        else:
            print 'Length.__mul__(): Cannot multiply a Length object by', other, \
                'Error return "False"'
            new = False

        return new

    def __div__(self, other):
        """Divides a Vector() object by a scalar.
            """
        if isinstance(other, int) or isinstance(other, float):
            new = Vector((array(self.get()) / other).tolist())
            new.units = self.units
        else:
            print 'Length.__div__(): Cannot divide a Length object by', other, \
                'Error return "False"'
            new = False

        return new


class datime():
    """Extend the datetime class to add mjd (Modified Julian Day) handling and
        some other conveniences.  I wanted to inherit from datetime, but that
        did not work because the datetime object, once created, is immutable,
        so I make a new date class with a datetime object dt as data, and mjd
        is calculated.
       
        Some examples: 
            d = datime()        sets d.dt to the current date and time
            d = datime(datetime(year, month[, day[, hour[, minute[, second[, us]
                                ]]]])) sets d.dt to given date/time 
            mjd = d.get()       returns the mjd (including fraction of a day)
            d.set(mjd)          will set the datetime in d.dt to the date 
                                corresponding to mjd
            d.set('2012-03-18 3:27:34.123','str')   will set the datetime 
                                in d.dt to the date in the string
            d.set('3:27:34.123','str')   will set the datetime in d.dt to 
                                today's date and the time in the string
            d.set('3:27','str') will set the datetime in d.dt to today's 
                                date and HH:MM in the string
            st = d.get('str')   sets st to a string YYYY-MM-DD HH:MM:SS.SSSSSS
            st = d.get('mstr')  sets st to a string YYYY-Mon-DD HH:MM:SS.SSSSSS
        
        Of course, all of the methods of datetime are available for date.dt also.
        """

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
            tz = 0  # Do not consider time zone correction.

        # Attach datetime() object to self
        self.dt = dtin
        # Calculate mjd day and add time as fraction of a day.  Zero of MJD is 1858/11/17.
        mjdt = self.dt - dt.datetime(1858, 11, 17, 0, 0, 0)
        mjd = mjdt.days + mjdt.seconds / 86400. + mjdt.microseconds / 86400000000.
        # Set our mjd (first call creates private attribute __mjd) 
        self.__mjd = mjd
        # Now correct for time zone so that time is in UTC.
        self.set(mjd + tz / 24.)

    def get(self, units='mjd'):
        """Return the date and time in various formats specified by units string
            'mjd' -- Modified Julian Day number, including fraction of a day (default)
            'str' -- standard datetime output string, YYYY-MM-DD HH:MM:SS.SSSSSS
            'mstr'  -- string with 3-char month text,  YYYY-Mon-DD HH:MM:SS.SSSSSS
            'trad'  -- time only, as an angle in radians"""

        if units == 'mjd':
            """Return the Modified Julian Day number, including fraction of a day"""
            return self.__mjd
        if units == 'str':
            """Output a standard Date/Time string (ms resolution, which can be truncated 
                by using slices if only date, or less time resolution is desired."""
            # Output string should have time to ms resolution, and the '.000' and [:23] is arranged to
            # guarantee this.
            return (self.dt.__str__() + '.000')[:23]
        if units == 'mstr':
            """Output standard Date/TIme string except with 3-char month abbreviation"""
            return (self.dt.strftime('%Y-%b-%d %H:%M:%S.%f') + '.000')[:24]
        elif units == 'tstamp':
            """Output LabVIEW timestamp (seconds since 1904-01-01)"""
            return (self.__mjd - 16480.0) * 86400.

    def set(self, value, units='mjd'):
        """Set the date and time, with input in various formats specified by units string
            'mjd' -- input is Modified Julian Day number, including fraction of a day (default)
            'str' -- standard datetime output string, YYYY-MM-DD HH:MM:SS.SSSSSS
            'mstr'  -- string with 3-char month text,  YYYY-Mon-DD HH:MM:SS.SSSSSS
            'trad'  -- time only, as an angle in radians"""
        if units == 'mjd':
            imjd = int(value)
            fmjd = (value - imjd) * 86400.
            self.dt = dt.datetime(1858, 11, 17, 0, 0, 0) + dt.timedelta(imjd, fmjd)
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
                us = int(ms) * 1000
            except:
                # Failure indicates that no fractional seconds were given,
                # so set microseconds to zero.
                ints = ss
                us = 0
            # Create a datetime object for this time, and initialize self with it.  Adding 1
            # us appears to be necessary to keep a time, e.g. 19:32:49.000000 from going to
            # 19:32:48.999999.
            self.__init__(dt.datetime(int(yr), int(mo), int(da), int(hh), int(mm), int(ints), us + 1))
        elif units == 'tstamp':
            mjd = value / 86400. + 16480.0  # Convert LabVIEW timestamp to mjd
            imjd = int(mjd)
            fmjd = (mjd - imjd) * 86400.
            self.dt = dt.datetime(1858, 11, 17, 0, 0, 0) + dt.timedelta(imjd, fmjd)
            self.__mjd = mjd
        return self.get()


from astropy.time import Time as astroTime
from astropy.time.core import TimeFromEpoch


class TimeLV(TimeFromEpoch):
    """
    LabVIEW timestamp input:
    number of seconds from 1904-01-01 00:00:00 UTC

    """
    # This corresponds to the zero reference time for LabVIEW.
    # Note that TAI and UTC are equivalent at the reference time.
    name = 'lv'
    unit = 1.0 / 86400.0  # Converts to seconds
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

        import pytz
        try:
            self.LocalTime = pytz.utc.localize(self.datetime, is_dst=None).astimezone(
                pytz.timezone('America/Los_Angeles'))
        except:
            try:
                locT = []
                for ll in self:
                    locT.append(
                        pytz.utc.localize(ll.datetime, is_dst=None).astimezone(pytz.timezone('America/Los_Angeles')))
                self.LocalTime = locT
            except:
                pass


from math import floor


class QuadraticInterpolator():
    """A quadratic interpolator class"""

    # Data
    def __init__(self, ytype='normal'):
        """Initialize the quadratic interpolator"""

        self.npt = 0  # Number of points currently in the interpolator (0, 1, 2 or 3)
        self.ytype = ytype  # Do not wrap values ['normal'] or wrap to -pi,pi ['signed'] or 0,2pi ['unsigned'])
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0  # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue)  # Initialize to empty

    # Methods
    def extendangle(self, a, b):
        """Given two angles A and B within the same 2.pi interval, return
           whichever of B - 2.pi, B, B + 2.pi is within pi of A."""

        # Convert to float, if possible [important if a, b are numpy arrays]
        a = float(a)
        b = float(b)
        dif = b - a
        if dif > pi:
            return b - 2 * pi
        elif dif < -pi:
            return b + 2 * pi
        else:
            return b

    def extend(self, x0, y0):
        """Pops the first value and pushes the given value (increments npts if needed), 
           then calls set()"""

        # Convert to float, if possible [important if x0, y0 are numpy arrays]
        x0 = float(x0)
        y0 = float(y0)
        # Remove the first value (superseded by the new value)
        self.x.pop(0)
        self.y.pop(0)

        self.npt += 1

        # If this is the first point entered, use it to initialize __X0
        if self.npt == 1:
            self.__X0 = x0
        elif self.npt > 3:
            # A valid point was popped off, so reset __X0 and x values for new minimum
            xdif = self.x[0]
            self.__X0 += xdif
            for i in range(0, 2):
                self.x[i] -= xdif

        self.x.append(x0 - self.__X0)  # Append the new x-coordinate after subtracting __X0
        self.y.append(y0)  # Append the new y-coordinate

        if self.npt > 3: self.npt = 3  # Make sure npt never exceeds 3

        self.set()  # Call the set() function to evaluate the coefficients

    def evaluate(self, x):
        """Given a time, evaluates the interpolation at that time"""

        twopi = 2 * pi
        c = self.c
        ytype = self.ytype

        # if ytype == 'normal': 
        x -= self.__X0

        y = x * (x * c[0] + c[1]) + c[2]

        if ytype == 'signed':
            return y - twopi * floor(y / twopi + 0.5)
        elif ytype == 'unsigned':
            return y - twopi * floor(y / twopi)
        else:
            return y

    def gradient(self, x):
        """Given a time, evaluates the gradient of the interpolation at that time"""

        # if self.ytype == 'normal':
        x -= self.__X0

        return 2.0 * x * self.c[0] + self.c[1]

    def empty(self, emptyvalue):
        """empty:   Clears the interpolator of all points, and sets npts to 0"""

        self.x = [0.0, 0.0, 0.0]  # The x coordinates (times) to use for the quadratic
        self.y = [0.0, 0.0, 0.0]  # The y coordinates corresponding to x
        self.c = [0.0, 0.0, 0.0]  # The coefficients calculated by set()
        self.npt = 0

    def canbracket(self, x):
        """canbracket: Says whether the interpolator can bracket the given time (true or false)"""

        # if self.ytype == 'normal': 
        x -= self.__X0

        if self.npt < 2:
            return False
        else:
            return (x >= self.x[3 - self.npt] and x <= self.x[self.npt - 1])

    def set(self):
        """Calculates the coefficients of the interpolator"""

        # Original routine does a sort of the x values at this point, and removes
        # duplicates, but this version assumes that x's are unique and in ascending
        # order, since if not, there is something else wrong.

        x = self.x
        y = self.y
        c = self.c

        if self.ytype != 'normal':
            if self.npt == 3:
                y[1] = self.extendangle(y[0], y[1])
                y[2] = self.extendangle(y[1], y[2])
            elif self.npt == 2:
                y[2] = self.extendangle(y[1], y[2])

        if self.npt == 0:
            # Return constant value self.emptyvalue
            c[0] = 0.0
            c[1] = 0.0
            c[2] = self.emptyvalue
        elif self.npt == 1:
            # Return constant value of y-coordinate
            c[0] = 0.0
            c[1] = 0.0
            c[2] = self.y[2]
        elif self.npt == 2:
            # Return coefficients for linear interpolation of two existing points
            c[0] = 0.0
            c[1] = (y[2] - y[1]) / (x[2] - x[1])
            c[2] = y[1] - c[1] * x[1]
        else:
            p = (y[2] - y[1]) / (x[2] - x[0]) / (x[2] - x[1])
            q = (y[1] - y[0]) / (x[2] - x[0]) / (x[1] - x[0])
            c[0] = p - q
            c[1] = q * (x[2] + x[1]) - p * (x[1] + x[0])
            c[2] = y[2] - c[0] * x[2] * x[2] - c[1] * x[2]

        # Copy the results
        for i in range(0, self.npt):
            self.c[i] = c[i]

    def getx(self):
        """Gets the contents of the interpolator"""

        x = []
        for i in range(len(self.x)):
            x.append(self.x[i] + self.__X0)
        return x


class RA_Interpolator(QuadraticInterpolator):
    """A convenience, to create an 'unsigned' quadratic interpolator as is needed for RA"""

    def __init__(self):
        """Initialize the RA quadratic interpolator"""
        self.npt = 0  # Number of points currently in the interpolator (0, 1, 2 or 3)
        self.ytype = 'unsigned'  # Wrap to 0,2pi ['unsigned']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0  # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue)  # Initialize to empty


class HA_Interpolator(QuadraticInterpolator):
    """A convenience, to create a 'signed' quadratic interpolator as is needed for 
       HA"""

    def __init__(self):
        """Initialize the HA quadratic interpolator"""
        self.npt = 0  # Number of points currently in the interpolator (0-3)
        self.ytype = 'signed'  # Wrap to -pi,pi ['signed']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0  # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue)  # Initialize to empty


class Dec_Interpolator(QuadraticInterpolator):
    """A convenience, to create a 'signed' quadratic interpolator as is needed for 
       Declination"""

    def __init__(self):
        """Initialize the quadratic interpolator"""
        self.npt = 0  # Number of points currently in the interpolator (0-3)
        self.ytype = 'signed'  # Wrap to -pi,pi ['signed']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0  # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue)  # Initialize to empty


import urllib2
from os import environ, path


def UT1_UTC(mjd):
    '''Download current IERS Bulletin A, locate and read UT1-UTC
       for given MJD day.  The time of day is ignored, which should
       lead to no more than 1.2 ms clock error (0.2 arcsec positional
       error, potentially corrected by phase calibration).
           Returns UT1-UTC as fraction of a day.
           Saves Bulletin as a file IERS.txt every 30 days, and falls
             back to saved file if URLError occurs
           Raises NameError if MJD cannot be found in Bulletin

       N.B.: Need to create an environment variable (not HOME) for
       placing the file.
       '''
    IERS_file = path.join(environ['HOME'], 'IERS.txt')
    #IERS_url = 'http://maia.usno.navy.mil/ser7/ser7.dat'
    IERS_url = 'https://datacenter.iers.org/data/latestVersion/6_BULLETIN_A_V2013_016.txt'

    try:
        h = urllib2.urlopen(IERS_url, timeout=1)
        # Copy the Bulletin IERS_file on disk, for potential fall-back
        with open(IERS_file, 'w') as f:
            f.write(h.read())
            h.close()
        # Safest to now open and read from the file, to make sure
        # the "fall-back" works.
        h = open(IERS_file, 'r')
    except:
        # DUT: Could not open URL for IERS Bulletin A, so read from
        #      existing file.
        print 'IERS Bulletin unreachable.  Reading from cached file.'
        h = open(IERS_file, 'r')

    for line in h:
        if line.find(str(int(mjd + 1))) == 19:
            # This actually uses "tomorrow's" value which is within 1 ms
            # of today, and anyway we tend to operate near the end of the
            # current day at OVRO.
            dut = float(line[49:62])
            h.close()
            break

    try:
        # Convert dut in seconds to fraction of a day
        dut = dut / 86400.
    except (UnboundLocalError, NameError):
        # This is likely to happen only if DUT is undefined, i.e. MJD
        # was not found in the file.  This should never happen, but...
        h.close()
        import os
        print "UT1_UTC: os.stat of IERS Bulletin A file\n",os.stat(IERS_file)
        print ""
        print "UT1_UTC: Today's MJD date not found in IERS Bulletin A."
        print "UT1_UTC: Setting DUT to zero"
        dut = 0

    return dut


def common_val_idx(array1, array2, precision=None):
    ''' Find the common values in two sorted arrays, and return the array
        of indexes of those common values in the two arrays.  If the
        parameter precision is given, the input arrays are rounded to
        that many decimal places before comparison.
        
    '''
    from numpy import intersect1d, searchsorted, round

    if precision:
        try:
            ar1 = round(array1 * 10 ** precision)
            ar2 = round(array2 * 10 ** precision)
        except Exception as e:
            print 'Error: common_val_idx: arrays could not be rounded to precision', precision
            print e
            return None, None
    else:
        ar1 = array1
        ar2 = array2
    common = intersect1d(ar1, ar2, True)
    idx1 = searchsorted(ar1, common)
    idx2 = searchsorted(ar2, common)
    return idx1, idx2


def nearest_val_idx(array1, array2):
    ''' Find the nearest values in the second array to the values in the first array
        and return the array of indexes of those nearest values in the second array.
    '''

    def find_nearest(array, value):
        from numpy import intersect1d, searchsorted
        from math import fabs
        idx = searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or fabs(value - array[idx - 1]) < fabs(value - array[idx])):
            return idx - 1
        else:
            return idx

    idx = []
    for value in array1:
        idx.append(find_nearest(array2, value))
    return array(idx)


def lobe(phi, mid=True):
    # Ensures that value phi lies between -pi and pi (if mid = True)
    # or 0 and 2*pi (if mid = False)
    intype = ''
    if isinstance(phi, list):
        intype = 'list'
        phi = array(phi)
    if isinstance(phi, float) or isinstance(phi, int):
        intype = 'float'
        phi = array([phi])
    if not isinstance(phi, ndarray):
        raise TypeError('Wrong argument type for lobe:' + str(type(phi)))
    phi = phi % (2 * pi)
    if mid:
        phi[where(phi > pi)] -= 2 * pi
        phi[where(phi < -pi)] += 2 * pi
    else:
        phi[where(phi < 0)] += 2 * pi
    if intype == 'list':
        return phi.tolist()
    elif intype == 'float':
        return phi[0]
    return phi


def ant_str2list(ant_str):
    ant_list = []
    try:
        grps = ant_str.split()
        for grp in grps:
            antrange = grp[3:].split('-')
            if len(antrange) == 1:
                if antrange != '':
                    ant_list.append(int(antrange[0]) - 1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0]) - 1, int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str', ant_str
        return None
    return array(ant_list)


def get_bl_order(n_ants):
    """Return the order of baseline data output by a CASPER correlator
    X engine."""
    order1, order2 = [], []
    for i in range(n_ants):
        for j in range(int(n_ants / 2), -1, -1):
            k = (i - j) % n_ants
            if i >= k:
                order1.append((k, i))
            else:
                order2.append((i, k))
    order2 = [o for o in order2 if o not in order1]
    return tuple([o for o in order1 + order2])


def bl_list(nant=16):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = nant*(nant-1)/2 + i.
    '''
    from numpy import ones
    bl2ord = ones((nant, nant), dtype='int') * (-1)
    nbl = nant * (nant - 1) / 2
    k = 0
    for i in range(nant - 1):
        for j in range(i + 1, nant):
            bl2ord[i, j] = k
            bl2ord[j, i] = k
            k += 1
    for i in range(nant):
        bl2ord[i, i] = nbl + i
    return bl2ord


bl2ord = bl_list()


def get_idbdir(t=None, usejsonfile=True):
    ''' Returns the root location of IDB files for the date given in Time object t.
        If t is not supplied, returns the root location for the latest data.
        
        Default in case of error is to return root location for the latest data.
    '''
    # Currently (2020-05-10), four variables are defined on pipeline:
    #   EOVSAJSON = /data1/eovsa/EOVSADB.json (if set, the other three are not needed)
    #   EOVSADB1 = /nas3/IDB/
    #   EOVSADB2 = /data1/eovsa/fits/IDB/
    #   EOVSADATE = <yyyy-mm-dd>  (date of start of data on EOVSADB2)
    import os
    from astropy.time import Time

    eovsajsonfile = os.getenv('EOVSADBJSON')
    if eovsajsonfile:
        if not os.path.exists(eovsajsonfile):
            print('GET_IDBDIR: JSON file {} does not exist. usejsonfile is set to False.'.format(eovsajsonfile))
            usejsonfile = False
    else:
        usejsonfile = False

    if usejsonfile:
        import json
        import numpy as np
        with open(eovsajsonfile) as f:
            eovsadb = json.load(f)
        dbarr = []
        datekey = eovsadb['EOVSADB'].keys()
        datetobj = Time(datekey)
        dtype = [('date', 'S10'), ('tobj', object), ('path', object)]
        for idx, k in enumerate(datekey):
            dbarr.append((k, datetobj[idx], eovsadb['EOVSADB'][k]))
        dbarr = np.array(dbarr, dtype=dtype)
        dbarrsort = np.sort(dbarr, order='tobj')

        if t is None:
            # Default to current time
            t = Time.now()
        print t.mjd, dbarrsort['tobj'][0].mjd
        if t.mjd < dbarrsort['tobj'][0].mjd:
            print('The date provided with t is before the time of the EOVSA first light. t is reset to {}.'.format(
                dbarrsort[0]['date']))
            t = dbarrsort['tobj'][0]
        # eodate = os.getenv('EOVSADATE')
        # Default to EOVSADB2
        # envar = 'EOVSADB2'
        datadir = dbarrsort['path'][np.where(t >= dbarrsort['tobj'])[0][-1]]
    else:
        if t is None:
            # Default to current time
            t = Time.now()
        eodate = os.getenv('EOVSADATE')
        # Default to EOVSADB2
        envar = 'EOVSADB2'
        if not eodate:
            # EOVSADATE is not defined
            print('GET_IDBDIR: Environment variable EOVSADATE is not defined. Returning root of latest data')
        else:
            try:
                if t < Time(eodate):
                    # Requested time is before EOVSADATE, so use EOVSADB1
                    envar = 'EOVSADB1'
            except:
                print('GET_IDBDIR: Invalid Time() object. Returning root of latest data.')
        datadir = os.getenv(envar)
    if not datadir:
        # Return default directory on pipeline
        datadir = '/data1/eovsa/fits/IDB/'
        print('GET_IDBDIR: Environment variable', envar, 'is not defined. Returning root of latest data.')
    if not datadir.endswith('/'):
        datadir = ''.join([datadir,'/'])
    return str(datadir)

def freq2bdname(fghz,t=None):
    '''Determine the band name from a given frequency in GHz, depending on date of observation.
       Just calls a different module.
    '''
    import chan_util_52 as cu52
    import chan_util_bc as cu34

    if t is None:
        t = Time.now()
    if t.mjd > 58536:
        return cu52.freq2bdname(fghz)
    else:
        return cu34.freq2bdname(fghz)

def fname2mjd(filename):
    ''' Get modified julian date from a standard IDB or UDB filename.
    
        Input can be a single string (filename) or a list or numpy array of filenames
        Output is a single mjd or a numpy array of mjds.
    '''
    from numpy import ndarray
    if type(filename) == ndarray or type(filename) == list:
        fstr = []
        for file in filename:
            fstem = file.split('/')[-1]
            fstr.append(fstem[3:7] + '-' + fstem[7:9] + '-' + fstem[9:11] + ' ' + fstem[11:13] + ':' + fstem[
                                                                                                       13:15] + ':' + fstem[
                                                                                                                      15:17])
        t = Time(fstr)
    else:
        fstem = filename.split('/')[-1]
        fstr = fstem[3:7] + '-' + fstem[7:9] + '-' + fstem[9:11] + ' ' + fstem[11:13] + ':' + fstem[
                                                                                              13:15] + ':' + fstem[
                                                                                                             15:17]
        t = Time(fstr)
    return t.mjd


def lin_phase_fit(f, pha, doplot=False):
    ''' Given an array of frequencies and corresponding phases,
        determine the best linear fit and return the parameters
        and standard deviation of the fit.  Optionally plots the
        phases and the fit, mainly for testing and evaluation
        purposes.  
        
        Inputs:
          f         Array of frequencies (does not need to be evenly spaced).
                       Note, no Nans allowed.
          pha       Array of phases, in radians, corresponding to array f. 
                       Note, any phases to be ignored can be flagged with Nan.
          doplot    Optional flag--if True, opens a new plot and plots the 
                       phases and the fit.
                       
        Returns:
          Numpy 3-element array of phase-offset, phase-slope, and 
          standard deviation of the fit
    '''
    import numpy as np
    if len(f) != len(pha):
        print 'Error: arrays not of same size:', len(f), len(pha)
        return None
    dpdf = []
    good, = np.where(np.logical_not(np.isnan(pha)))
    if len(good) < 3:
        # Not enough points to fit, so return zeros and a large standard deviation
        return np.array((0, 0, np.pi))
    f_ = f[good]
    pha_ = pha[good]
    for i in range(len(f_) - 1):
        dpdf.append((pha_[i + 1] - pha_[i]) / (f_[i + 1] - f_[i]))
    dpdf = np.array(dpdf)
    slp = np.median(dpdf)
    p = np.polyfit(f_, np.unwrap((pha_ - f_ * slp)), 1)
    p[0] = p[0] + slp
    stdev = np.std(lobe(pha_ - np.polyval(p, f_)))
    if doplot:
        import matplotlib.pylab as plt
        plt.plot(f, pha, '.')
        plt.plot(f, lobe(np.polyval(p, f)))
    return np.array((p[1], p[0], stdev))


def fix_time_drift(out):
    ''' Routine to correct a linear phase drift vs. time for baselines with
        Ant 14 in a standard read_idb() output file.  This calculates a slope vs.
        time for each frequency and pol = XX and YY separately, but then uses 
        the median of the delay (slope/fghz) for well-determined slopes (those 
        with stdev < 0.7) to correct the phase on ALL frequencies and polarizations. 
        Other baselines than those with Ant 14 are returned unmodified.
        (Although they COULD be corrected -- TODO)
    '''
    import numpy as np
    nant, npol, nband, nt = out['x'][bl2ord[13, :13]].shape  # Consider baselines with Ant14 only
    for iant in range(nant):
        for ipol in range(2):  # Use only polarizations XX and YY for slope determination
            slopes = []
            for iband in range(nband):
                phz = np.angle(out['x'][bl2ord[13, iant], ipol, iband])
                # if out['flags'][iant,ipol,iband] == 0:
                p = lin_phase_fit(out['time'], phz)
                if p[2] < 0.7:
                    slopes.append(p[1] / out['fghz'][iband])

        if len(slopes) > 0:
            dpdt = np.nanmedian(slopes)  # Radians/GHz/Day
            for ipol in range(npol):  # Apply the corrections to all polarization products
                for iband in range(nband):
                    pfit = dpdt * out['fghz'][iband] * (out['time'] - out['time'][nt / 2])
                    out['x'][bl2ord[13, iant], ipol, iband] *= np.cos(pfit) - 1j * np.sin(pfit)
    return out


def sat_phase(out, ant, doplot=False):
    ''' Takes the autocorrelation read from a packet dump (PRT) file taken while the array is
        observing an R/L polarized communications satellite (e.g. Ciel), for one time sample, and
        corrects the XY and YX phase for phase slope (X-Y delay) and offset.  The correction
        is done "in place," so nothing is returned, but the "out" array is corrected.  This
        works only for one antenna at a time (specified by "ant"),  If doplot is True, it
        also plots the corrected phase (mainly for debugging purposes).

        Inputs:    
        out     is out['a'][:,:,:,n], from pcapture2.rd_jspec(file), where n is the index of
                   one time sample, and file is the name of a PRT (packet dump) file.
        ant     index of the antenna for which to do the correction.           
        
        Optional parameters:
        doplot  if True, plots the corrected phase vs. channel
        
        Ultimately, this routine should return the phase slope and offset, which could form
        the basis of a calibration scheme, but this is still under investigation. It should
        also correct and return the X vs. Y amplitude scaling.
    '''
    import numpy as np
    xy = out[ant, 2, :]
    yx = out[ant, 3, :]
    f = np.arange(len(out[0, 0]))
    phi = lobe(np.angle(xy) - np.angle(yx))
    good = np.where(np.abs(xy) > np.median(np.abs(xy)))
    poff, pslp, pstd = lin_phase_fit(f[good], phi[good])
    if np.mean(lobe(np.angle(xy[1100:1300]) - (pslp * f[1100:1300] + poff + pi) / 2.)) > 0: poff = poff + 2 * np.pi
    if doplot:
        import matplotlib.pylab as plt
        plt.plot(f, lobe(np.angle(xy) - (pslp * f + poff + np.pi) / 2.), '.')
        plt.plot(f, lobe(np.angle(yx) + (pslp * f + poff + np.pi) / 2.), '.')
    out[ant, 2, :] *= np.exp(-1j * (pslp * f + poff + np.pi) / 2.)
    out[ant, 3, :] *= np.exp(1j * (pslp * f + poff + np.pi) / 2.)


def plot_sched_log(tail=None):
    ''' Routine to read the schedule.log file on Helios and plot the lines giving
        duration anomalies.  Normally inc_time in the schedule should be run every
        second with a jitter of less than 10 ms.  Lines with non-standard delays
        are written to the log file.  This routine simply plots that information,
        to get an idea of when something is causing delays.
        
        keywords:
        
        tail    A number N that specifies that only the last N (or fewer) durations 
                   should be plotted.  If omitted, all duration lines are plotted.
        
        Note: This routine MUST be run on Helios (since that is where the log file is).
    '''
    import matplotlib.pylab as plt
    from util import Time
    f = open('/tmp/schedule.log', 'r')
    lines = f.readlines()
    f.close()
    pd_list = []
    dur_list = []
    for line in lines:
        tok = line.strip().split()
        if len(tok) == 3:
            try:
                t = Time(tok[0] + ' ' + tok[1]).plot_date
                pd_list.append(t)
                dur_list.append(float(tok[2]))
            except:
                pass
    if tail is None:
        plt.plot_date(pd_list, dur_list, '.')
    else:
        plt.plot_date(pd_list[-tail:], dur_list[-tail:], '.')
