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
# * 

import StringUtil as su
from numpy import pi, sqrt, array, mat, matrix, dot, where, ndarray
import datetime as dt
from time import gmtime

class Angle:
    """General angle class, converts input to radians, but handles input
       and output in various units."""

    __pi               = pi
    __twopi            = 2*pi
    __degperrad        = 180./pi
    __arcsecperdeg     = 3600.
    __masperdeg        = 3600000.
    __arcsecperrad     = __arcsecperdeg * __degperrad
    __arcminperrad     = 60 * __degperrad
    __masperrad        = __masperdeg * __degperrad
    
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
                strg = ' 0'+strg
            elif len(strg) == 12:
                # Need to add missing zero or space
                if strg[0] == '-':
                    strg = strg[0]+'0'+strg[1:]
                else:
                    strg = ' '+strg
            return strg
        elif units == 'hmsstr' or units == 'hms':
            # String in the form -HH:MM:SS.SSS
            strg = su.dmsStrFromDeg(self.radians * self.__degperrad / 15., precision=3)
            if len(strg) == 11:
                # Need to add missing zero
                strg = '0'+strg
            return strg
        elif units == 'masec':
            return self.radians * self.__masperrad
        elif units == 'arcsec':
            return self.radians * self.__arcsecperrad
        else:
            print 'ERROR: Unknown units string, radians used'
            return self.radians

    def __add__(self,other):
        """Add one Angle object to another, or add a scalar to
            an Angle object.
            """
        if isinstance(other,int) or isinstance(other, float):
            # scalar is interpreted as radians
            radians = other
        elif isinstance(other,Angle):
            radians = other.radians
        else:
            print 'Angle.__add__(): Cannot add', other, 'to an Angle()'
            return 0
        
        # Create the right kind of object to return (same as self)
        if isinstance(self,RA_Angle):
            return RA_Angle(self.radians + radians)
        elif isinstance(self, Dec_Angle):
            return Dec_Angle(self.radians + radians)
        else:
            return Angle(self.radians + radians)

    def __sub__(self,other):
        """Subtract one Angle object from another, or subtract a scalar from
            an Angle object.
            """
        if isinstance(other,int) or isinstance(other, float):
            # scalar is interpreted as radians
            radians = other
        elif isinstance(other,Angle):
            radians = other.radians
        else:
            print 'Angle.__sub__(): Cannot subtract', other, 'from an Angle()'
            return 0
        
        # Create the right kind of object to return (same as self)
        if isinstance(self,RA_Angle):
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
            twopi = 2.*pi
            while radians >= pi:
                radians -= twopi
            while radians < -pi:
                radians += twopi
            self.radians = radians
            # Set the decvalid flag to False if not between -pi/2 and pi/2.
            if radians >= pi/2 or radians < -pi/2:
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
            print 'Length.set(): Invalid length value',value,'Setting to zero'
            self.meters = 0.0
            
        if units == 'cm':
            self.meters = value/cmperm
        elif units == 'm':
            self.meters = float(value)
        elif units == 'km':
            self.meters = value*mperkm
        else:
            print 'Length.set(): Invalid length units',units,'Meters assumed.'

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
            value = self.meters*cmperm
        elif units == 'm':
            value = self.meters
        elif units == 'km':
            value = self.meters/mperkm
            
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
        if isinstance(other,Length):
            new = Length(self.meters + other.meters)
        elif isinstance(other,int) or isinstance(other,float):
            new = Length(self.meters + other)
        else:
            print 'Length.__add__(): Cannot add',other,\
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
        if isinstance(other,Length):
            new = Length(self.meters - other.meters)
        elif isinstance(other,int) or isinstance(other,float):
            new = Length(self.meters - other)
        else:
            print 'Length.__sub__(): Cannot subtract',other,\
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
        if isinstance(other,Length):
            new = self.meters / other.meters
        elif isinstance(other,int) or isinstance(other,float):
            new = Length(self.meters / other)
            new.units = self.units
        else:
            print 'Length.__div__(): Cannot divide a Length object by',other,\
                'Error return "False"'
            new = False
        
        return new

class Vector():
    """A general class to hold and manipulate 3-element vectors of
        Length() objects
        """

    def __init__(self, lengths = None, units='m'):
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
        self.set(lengths,units)

    def set(self, lengths = None, units='m'):
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
            if not isinstance(lengths[i],Length):
                # One of the supplied objects is not a Length() object, so
                # see if it is a scalar constant.
                if isinstance(lengths[i],int) or isinstance(lengths[i],float):
                    # It is, so create a Length() object using it and the
                    # supplied units.
                    lengths[i] = Length(lengths[i],units=units)
                else:
                    # It is neither a Length object nor a scalar constant, so
                    # it is an error.
                    print 'Vector.set(): Item',i,'is not a Length() or a constant.\n\
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
            #if units == 'str':
            #    value = self.lengths[i].get(units)
            #    out.append(str(value) + ' ' + self.units)
            #else:
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
            print 'Vector.__add__(): Object',other,'must be a Vector() object.',\
                'Returning a zero vector.'
            return Vector([0.0,0.0,0.0])
        
    def __sub__(self, other):
        """Subtract two vector objects"""
        if isinstance(other, Vector):
            out = Vector()
            for i in range(3):
                out.lengths[i] = self.lengths[i] - other.lengths[i]
            return out
        else:
            print 'Vector.__sub__(): Object',other,'must be a Vector() object.',\
                'Returning a zero vector.'
            return Vector([0.0,0.0,0.0])
        
    def magnitude(self, units='m'):
        """Returns the magnitude of the Vector in the units specified.
            """
        if units == 'str':
            # Return the lengths in the units of the Vector object
            vals = self.get(self.units)
        else:
            # Return the lengths in the specified units
            vals = self.get(units)
            
        #Calculate the magnitude
        mag = 0.0
        for i in range(3):
            mag += vals[i]*vals[i]
        magnitude = sqrt(mag)
        if units == 'str':
            # Make a string out of it
            return str(magnitude)+' '+self.units
        
        return magnitude

    def rotate(self, rot=None):
        """Rotate a Vector object by applying the rotation matrix rot."""

        if rot is None:
            rot = mat([[1,0,0],[0,1,0],[0,0,1]])
        if not isinstance(rot,matrix):
            print 'Vector.rotate(): Object',rot,'is not a matrix.  No rotation applied.'
            rot = mat([[1,0,0],[0,1,0],[0,0,1]])

        # Apply rotation (output is in meter)
        vnew = dot(rot,self.get())
        # Create a new Vector object for output.  Have to convert vnew numpy
        # array to a "list of lists", then take the first element
        new = Vector(vnew.tolist()[0])
        new.units = self.units  # Set output units same as input
        return new

    def __mul__(self, other):
        """Multiplies a Vector() object by a scalar.
            """
        if isinstance(other,int) or isinstance(other,float):
            new = Vector((array(self.get()) * other).tolist())
            new.units = self.units
        else:
            print 'Length.__mul__(): Cannot multiply a Length object by',other,\
                'Error return "False"'
            new = False
        
        return new
    def __div__(self, other):
        """Divides a Vector() object by a scalar.
            """
        if isinstance(other,int) or isinstance(other,float):
            new = Vector((array(self.get()) / other).tolist())
            new.units = self.units
        else:
            print 'Length.__div__(): Cannot divide a Length object by',other,\
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
            tz = 0   # Do not consider time zone correction.

        # Attach datetime() object to self
        self.dt = dtin
        # Calculate mjd day and add time as fraction of a day.  Zero of MJD is 1858/11/17.
        mjdt = self.dt - dt.datetime(1858,11,17,0,0,0)
        mjd = mjdt.days + mjdt.seconds/86400. + mjdt.microseconds/86400000000.
        # Set our mjd (first call creates private attribute __mjd) 
        self.__mjd = mjd
        # Now correct for time zone so that time is in UTC.
        self.set(mjd + tz/24.)

    def get(self,units='mjd'):
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
            return (self.__mjd - 16480.0)*86400.

    def set(self, value, units='mjd'):
        """Set the date and time, with input in various formats specified by units string
            'mjd' -- input is Modified Julian Day number, including fraction of a day (default)
            'str' -- standard datetime output string, YYYY-MM-DD HH:MM:SS.SSSSSS
            'mstr'  -- string with 3-char month text,  YYYY-Mon-DD HH:MM:SS.SSSSSS
            'trad'  -- time only, as an angle in radians"""
        if units == 'mjd':
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

        import pytz
        try:
            self.LocalTime = pytz.utc.localize(self.datetime, is_dst=None).astimezone(pytz.timezone('America/Los_Angeles'))
        except:
            try:
                locT = []
                for ll in self:
                    locT.append(pytz.utc.localize(ll.datetime, is_dst=None).astimezone(pytz.timezone('America/Los_Angeles')))
                self.LocalTime = locT
            except:
                pass


from math import floor

class QuadraticInterpolator():

    """A quadratic interpolator class"""

    # Data
    def __init__(self, ytype='normal'):
        """Initialize the quadratic interpolator"""
        
        self.npt = 0           # Number of points currently in the interpolator (0, 1, 2 or 3)
        self.ytype = ytype  # Do not wrap values ['normal'] or wrap to -pi,pi ['signed'] or 0,2pi ['unsigned'])
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0        # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue) # Initialize to empty

    # Methods
    def extendangle(self, a, b):
        """Given two angles A and B within the same 2.pi interval, return
           whichever of B - 2.pi, B, B + 2.pi is within pi of A."""

        # Convert to float, if possible [important if a, b are numpy arrays]
        a = float(a)
        b = float(b)
        dif = b - a
        if dif > pi:
            return b - 2*pi
        elif dif < -pi:
            return b + 2*pi
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
        if self.npt == 1: self.__X0 = x0
        elif self.npt > 3:
            # A valid point was popped off, so reset __X0 and x values for new minimum
            xdif = self.x[0]
            self.__X0 += xdif
            for i in range(0,2):
                self.x[i] -= xdif

        self.x.append(x0-self.__X0)   # Append the new x-coordinate after subtracting __X0
        self.y.append(y0)             # Append the new y-coordinate

        if self.npt > 3: self.npt = 3  # Make sure npt never exceeds 3

        self.set()    # Call the set() function to evaluate the coefficients

    def evaluate(self,x):
        """Given a time, evaluates the interpolation at that time"""

        twopi = 2*pi
        c = self.c
        ytype = self.ytype

        # if ytype == 'normal': 
        x -= self.__X0

        y = x * (x * c[0] + c[1]) + c[2]

        if ytype == 'signed': 
            return y - twopi * floor(y/twopi + 0.5)
        elif ytype == 'unsigned':
            return y - twopi * floor(y/twopi)
        else:
            return y
        
    def gradient(self,x): 
        """Given a time, evaluates the gradient of the interpolation at that time"""

        #if self.ytype == 'normal': 
        x -= self.__X0

        return 2.0 * x * self.c[0] + self.c[1]

    def empty(self,emptyvalue):
        """empty:   Clears the interpolator of all points, and sets npts to 0"""

        self.x = [0.0, 0.0, 0.0]   # The x coordinates (times) to use for the quadratic
        self.y = [0.0, 0.0, 0.0]   # The y coordinates corresponding to x
        self.c = [0.0, 0.0, 0.0]   # The coefficients calculated by set()
        self.npt = 0

    def canbracket(self,x):
        """canbracket: Says whether the interpolator can bracket the given time (true or false)"""

        # if self.ytype == 'normal': 
        x -= self.__X0

        if self.npt < 2: 
            return False
        else:
            return (x >= self.x[3-self.npt] and x <= self.x[self.npt-1])

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
            c[2] = y[1] - c[1]*x[1]
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
            x.append(self.x[i]+self.__X0)
        return x

class RA_Interpolator(QuadraticInterpolator):
    """A convenience, to create an 'unsigned' quadratic interpolator as is needed for RA"""
    def __init__(self):
        """Initialize the RA quadratic interpolator"""
        self.npt = 0           # Number of points currently in the interpolator (0, 1, 2 or 3)
        self.ytype = 'unsigned'  # Wrap to 0,2pi ['unsigned']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0        # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue) # Initialize to empty

class HA_Interpolator(QuadraticInterpolator):
    """A convenience, to create a 'signed' quadratic interpolator as is needed for 
       HA"""
    def __init__(self):
        """Initialize the HA quadratic interpolator"""
        self.npt = 0           # Number of points currently in the interpolator (0-3)
        self.ytype = 'signed'  # Wrap to -pi,pi ['signed']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0        # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue) # Initialize to empty

class Dec_Interpolator(QuadraticInterpolator):
    """A convenience, to create a 'signed' quadratic interpolator as is needed for 
       Declination"""
    def __init__(self):
        """Initialize the quadratic interpolator"""
        self.npt = 0           # Number of points currently in the interpolator (0-3)
        self.ytype = 'signed'  # Wrap to -pi,pi ['signed']
        self.emptyvalue = 0.0  # The value to use if npts = 0
        self.__X0 = 0.0        # Value to be subtracted from each x value to maintain accuracy
        self.empty(self.emptyvalue) # Initialize to empty


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
    IERS_file = path.join(environ['HOME'],'IERS.txt')
    IERS_url = 'http://maia.usno.navy.mil/ser7/ser7.dat'


    try:
        h = urllib2.urlopen(IERS_url,timeout=1)
        # Copy the Bulletin IERS_file on disk, for potential fall-back
        with open(IERS_file,'w') as f:
            f.write(h.read())
            h.close()
        # Safest to now open and read from the file, to make sure
        # the "fall-back" works.
        h = open(IERS_file,'r')                
    except:
        # DUT: Could not open URL for IERS Bulletin A, so read from
        #      existing file.
        print 'IERS Bulletin unreachable.  Reading from cached file.'
        h = open(IERS_file,'r')

    for line in h:
        if line.find(str(int(mjd+1))) == 19:
            # This actually uses "tomorrow's" value which is within 1 ms
            # of today, and anyway we tend to operate near the end of the
            # current day at OVRO.
            dut = float(line[49:62])
            h.close()
            break

    try:
        # Convert dut in seconds to fraction of a day
        dut = dut/86400.
    except (UnboundLocalError, NameError):
        # This is likely to happen only if DUT is undefined, i.e. MJD
        # was not found in the file.  This should never happen, but...
        h.close()
        raise NameError("UT1_UTC: Today's MJD date not found in IERS Bulletin A.")

    return dut

def common_val_idx(array1,array2):
    ''' Find the common values in two sorted arrays, and return the array
        of indexes of those common values in the two arrays.
    '''
    from numpy import intersect1d,searchsorted
    
    common = intersect1d(array1,array2,True)
    idx1 = searchsorted(array1,common)
    idx2 = searchsorted(array2,common)
    return idx1, idx2

def nearest_val_idx(array1,array2):
    ''' Find the nearest values in the second array to the values in the first array
        and return the array of indexes of those nearest values in the second array.
    '''
    def find_nearest(array,value):
        from numpy import intersect1d,searchsorted
        from math import fabs
        idx = searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or fabs(value - array[idx-1]) < fabs(value - array[idx])):
            return idx-1
        else:
            return idx
    idx = []
    for value in array1:
        idx.append(find_nearest(array2,value))
    return array(idx)

def lobe(phi,mid=True):
    # Ensures that value phi lies between -pi and pi (if mid = True)
    # or 0 and 2*pi (if mid = False)
    intype = ''
    if isinstance(phi,list):
        intype = 'list'
        phi = array(phi)
    if isinstance(phi,float) or isinstance(phi,int):
        intype = 'float'
        phi = array([phi])
    if not isinstance(phi,ndarray):
        raise TypeError('Wrong argument type for lobe:'+str(type(phi)))
    phi = phi % (2*pi)
    if mid:
        phi[where(phi > pi)] -= 2*pi
        phi[where(phi < -pi)] += 2*pi
    else:
        phi[where(phi < 0)] += 2*pi
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
                    ant_list.append(int(antrange[0])-1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0])-1,int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str',ant_str
        return None
    return array(ant_list)
