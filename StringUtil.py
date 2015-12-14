#!/usr/bin/env python
"""String utilities, with an emphasis on support for sexagesimal numbers
(e.g. degrees:minutes:seconds).

To do:
- Add tests for dmsStrFromDeg with omitExtraFields true,
  dmsStrFromSec, quoteStr etc.

History:
2001-03-12 ROwen    strToFloat modified to handle ".".
2001-03-14 ROwen    Added prettyDict; still need to add a test for it.
2001-04-13 ROwen    Fixed a bug in splitDMSStr: -0:... was not retaining the minus sign;
    added dmsStrToSec, which may be useful for handling durations;
    modified splitDMSStr to allow an arbitrary number digits in the first field
    (use range checking of the equivalent numeric value to limit this),
    to return sing as a separate entity, and to return all strings;
    improved the test suite, including adding a nearly silent test.
2001-04-20 ROwen    Fixed a bug in degToDMSStr: borderline cases could print
    as 60 instead of 00 and incrementing the next field.
2001-07-22 ROwen    Made DegStr work on unix as well as Mac
    (I still don't know the Windows character).
2001-08-29 ROwen    Made DegStr and DMSStr unicode strings, a reasonable choice
    now that one can print unicode characters.
2001-11-01 ROwen    Added nFields argument to degToDMSStr and secToDMSStr;
    changed degToDMSStr, secToDMSStr and neatenDMSStr to eliminate leading spaces.
2002-08-08 ROwen    Moved to RO and renamed from ROStringUtil. Renamed
     xToY functions to yFromX. Hid private constants with leading underscore.
2002-08-21 ROwen    Modified splitDMSStr to require all but the first field be <60.
2003-04-04 ROwen    Bug fix: dmsStrFromDeg and dmsStrFromSec mis-handled
    some negative values.
2003-06-19 ROwen    Fixed two tests which broke in upgrades listed above.
2003-07-15 ROwen    Added omitExtraFields to degFromDMSStr and secFromDMSStr;
                    added quoteStr function.
2003-10-27 ROwen    Added plural function.
2004-01-09 ROwen    Added AngstromStr, LambdaStr and MuStr constants.
2004-05-18 ROwen    Stopped importing sys since it was not being used.
2005-06-27 ROwen    Fixed a nonfunctional assert statement in splitDMSStr.
2006-01-09 ROwen    Fixed a bug in dmsStrFromDeg: dmsStrFromDeg(-50.650002) = "-50:38:60.0"
                    and a related bug in dmsStrFroMSec. Thanks to John Lucey.
                    (Note: the test code had a case for it, but expected the wrong value.)
2007-06-04 ROwen    Bug fix: dmsStrFromSec gave bad results if nFields != 3.
2008-01-30 ROwen    Removed unused variable signNum (found by pychecker).
2008-04-29 ROwen    Added strFromException, a unicode-safe replacement for str(exception).
2008-05-02 ROwen    Made prettyDict unicode-safe by using repr.
2008-11-14 ROwen    Added unquoteStr.
"""
import re

AngstromStr = u"\N{ANGSTROM SIGN}"
DegStr = u"\N{DEGREE SIGN}"
DMSStr = DegStr + u"'\""
LambdaStr = u"\u00c5" # for some reason this fails: u"\N{GREEK SMALL LETTER LAMBDA}"
MuStr = u"\N{GREEK SMALL LETTER MU}"

def dmsStrFromDeg (decDeg, nFields=3, precision=1, omitExtraFields = False):
    """Convert a number to a sexagesimal string with 1-3 fields.

    Inputs:
    - decDeg: value in decimal degrees or hours
    - nFields: number of fields; <=1 for dddd.ddd, 2 for dddd:mm.mmm, >=3 for dddd:mm:ss.sss
    - precision: number of digits after the decimal point in the last field;
        if 0, no decimal point is printed; must be >= 0
    - omitExtraFields: omit fields that are zero, starting from the right
    
    Error conditions:
    - Raises ValueError if precision < 0
    """
    nFields = min(3, nFields)
    signStr, fieldStrs = _getDMSFields(decDeg, nFields, precision)
    
    if omitExtraFields:
        while fieldStrs and float(fieldStrs[-1]) == 0.0:
            fieldStrs.pop(-1)
            
    return signStr + ":".join(fieldStrs)

def dmsStrFromSec (decSec, nFields=3, precision=1, omitExtraFields = True):
    """Convert a number, in seconds, to a sexagesimal string.
    Similar to dmsStrFromDeg, but takes seconds, not degrees,
    and omitExtraFields omits fields from the left, not the right.

    Inputs:
    - decSec: value in decimal seconds or arc seconds
    - nFields: number of fields; <=1 for ss.sss, 2 for mm:ss.ss, >= 3 for dddd:mm:ss.sss
    - precision: number of digits after the decimal point in the seconds field;
        if 0, no decimal point is printed; must be >= 0
    - omitExtraFields: omit fields that are zero, starting from the left.
    
    Error conditions:
    - Raises ValueError if precision < 0
    """
    nFields = min(3, nFields)
    if nFields < 1:
        raise ValueError("nFields=%r; must be >= 1" % (nFields,))
    adjNum = decSec / (60.0**(nFields-1))
    signStr, fieldStrs = _getDMSFields(adjNum, nFields, precision)
    
    if omitExtraFields:
        while fieldStrs and float(fieldStrs[0]) == 0.0:
            fieldStrs.pop(0)
            
    return signStr + ":".join(fieldStrs)

def degFromDMSStr (dmsStr):
    """Convert a string of the basic form dd[:mm[:ss]] to decimal degrees.
    See splitDMSStr for details of the format.
    
    Error conditions:
    - Raises ValueError if the string cannot be parsed
    """
    dmsItems = splitDMSStr(dmsStr)

    # extract sign
    if dmsItems[0] == '-':
        signMult = -1.0
    else:
        signMult = 1.0
    dmsItems[0:1] = []

    # combine last two elements and convert to float
    dmsItems[-2:] = [floatFromStr(dmsItems[-2]) + floatFromStr(dmsItems[-1])]

    # convert all but last item to float
    for ind in range(len(dmsItems) - 1):
        dmsItems[ind] = intFromStr(dmsItems[ind])
    
    dmsItems.reverse()
    decDeg = 0.0
    for dmsField in dmsItems:
        decDeg = abs(dmsField) + (decDeg / 60.0)
    return signMult * decDeg

def secFromDMSStr (dmsStr):
    """Convert a string of the basic form [[dd:]mm:]ss to decimal degrees.
    Note that missing fields are handled differently than degFromDMSStr!
    See splitDMSStr for details of the format.
    
    error conditions:
        raises ValueError if the string cannot be parsed
    """
    dmsItems = splitDMSStr(dmsStr)

    # extract sign
    if dmsItems[0] == '-':
        signMult = -1.0
    else:
        signMult = 1.0
    dmsItems[0:1] = []

    # combine last two elements and convert to float
    dmsItems[-2:] = [floatFromStr(dmsItems[-2]) + floatFromStr(dmsItems[-1])]

    # convert all but last item to float
    for ind in range(len(dmsItems) - 1):
        dmsItems[ind] = intFromStr(dmsItems[ind])
    
    decSec = 0.0
    for dmsField in dmsItems:
        decSec = abs(dmsField) + (decSec * 60.0)
    return signMult * decSec

def secStrFromDMSStr(dmsStr):
    """Convert a string of the basic form [[dd:]mm:]ss to decimal seconds
    preserving the original accuracy of seconds
    Note that missing fields are handled differently than degFromDMSStr!
    See splitDMSStr for details of the format.
    
    error conditions:
        raises ValueError if the string cannot be parsed
    """
    dmsItems = splitDMSStr(dmsStr)

    # extract sign and fractional seconds (includes decimal point)
    signStr = dmsItems[0]
    fracSecStr = dmsItems[-1]

    # compute integer seconds
    # convert all but first and last items to integers
    intList = [intFromStr(item) for item in dmsItems[1:-1]]
    
    intSec = 0
    for intVal in intList:
        intSec = abs(intVal) + (intSec * 60)
    return "%s%s%s" % (signStr, intSec, fracSecStr)

FloatChars = "0123456789+-.eE"

def checkDMSStr(dmsStr):
    """Verify a sexagesimal string; returns True if valid, False if not
    """
    try:
        splitDMSStr(dmsStr)
        return True
    except:
        return False

def dmsStrFieldsPrec(dmsStr):
    """Return the following information about a sexagesimal string:
    - the number of colon-separated fields
    - the precision of the right-most field
    """
    if dmsStr == "":
        return (0, 0)

    precArry = dmsStr.split(".")
    if len(precArry) > 1:
        precision = len(precArry[1])
    else:
        precision = 0
    nFields = dmsStr.count(":") + 1
    return (nFields, precision)
    
def findLeftNumber(astr, ind):
    """Find the starting and ending index of the number
    enclosing or to the left of index "ind".
    Return (None, None) if no number found.

    Warning: this is not a sophisticated routine. It looks for
    the a run of characters that could be present
    in a floating point number. It does not sanity checking
    to see if they make a valid number.
    """
    leftInd = _findLeftOfLeftNumber(astr, ind)
    if leftInd == None:
        return (None, None)
    rightInd = _findRightOfRightNumber(astr, leftInd)
    if rightInd == None:
        return (None, None)
    return (leftInd, rightInd)

def findRightNumber(astr, ind):
    """Find the starting and ending index of the number
    enclosing or to the right of index "ind".
    Returns (None, None) if no number found.

    Warning: this is not a sophisticated routine. It looks for
    the a run of characters that could be present
    in a floating point number. It does not sanity checking
    to see if they make a valid number.
    """
    rightInd = _findRightOfRightNumber(astr, ind)
    if rightInd == None:
        return (None, None)
    leftInd = _findLeftOfLeftNumber(astr, rightInd)
    if leftInd == None:
        return (None, None)
    return (leftInd, rightInd)

def _findLeftOfLeftNumber(astr, ind):
    """Find the index of the first character of the number
    enclosing or to the left of index "ind".
    Returns None if no number found.

    Warning: this is not a sophisticated routine. It looks for
    the left-most of a run of characters that could be present
    in a floating point number. It does not sanity checking
    to see if they make a valid number.
    """
    leftInd = None
    for tryind in range(ind, -1, -1):
        if astr[tryind] in FloatChars:
            leftInd = tryind
        elif leftInd != None:
            break
    return leftInd

def _findRightOfRightNumber(astr, ind):
    """Find the index of the last character of the number
    enclosing or to the right of index "ind".
    Returns None if no number found.

    Warning: this is not a sophisticated routine. It looks for
    the right-most of a run of characters that could be present
    in a floating point number. It does not sanity checking
    to see if they make a valid number.
    """
    rightInd = None
    for tryind in range(ind, len(astr)):
        if astr[tryind] in FloatChars:
            rightInd = tryind
        elif rightInd != None:
            break
    return rightInd

def neatenDMSStr (dmsStr):
    """Convert a sexagesimal string to a neater version.
    
    error conditions:
        raises ValueError if the string cannot be parsed
    """
    if dmsStr == "":
        return ""

    precArry = dmsStr.split(".")
    if len(precArry) > 1:
        precision = len(precArry[1])
    else:
        precision = 0
    fieldArry = dmsStr.split(":")
    nFields = len(fieldArry)
    
    floatValue = degFromDMSStr(dmsStr)
    return dmsStrFromDeg(floatValue, nFields=nFields, precision=precision)

def plural(num, singStr, plStr):
    """Return singStr or plStr depending if num == 1 or not.
    A minor convenience for formatting messages (in lieu of ?: notation)
    """
    if num == 1:
        return singStr
    return plStr

def prettyDict(aDict, entrySepStr = "\n", keyValSepStr = ": "):
    """Format a dictionary in a nice way
    
    Inputs:
    aDict: the dictionary to pretty-print
    entrySepStr: string separating each dictionary entry
    keyValSepStr: string separating key and value for each entry
    
    Returns a string containing the pretty-printed dictionary
    """
    sortedKeys = aDict.keys()
    sortedKeys.sort()
    eltList = []
    for aKey in sortedKeys:
        eltList.append(repr(aKey) + keyValSepStr + repr(aDict[aKey]))
    return entrySepStr.join(eltList)

# constants used by splitDMSStr
# DMSRE = re.compile(r"^\s*([+-]?)(\d{0,3})\s*(?:\:\s*(\d{0,2})\s*){0,2}(\.\d*)?\s*$")
_DegRE =       re.compile(r"^\s*([+-]?)(\d*)(\.\d*)?\s*$")
_DegMinRE =    re.compile(r"^\s*([+-]?)(\d*)\s*\:\s*([0-5]?\d?)(\.\d*)?\s*$")
_DegMinSecRE = re.compile(r"^\s*([+-]?)(\d*)\s*\:\s*([0-5]?\d?):\s*([0-5]?\d?)(\.\d*)?\s*$")

def splitDMSStr (dmsStr):
    """Split a sexagesimal string into fields
    returns one of the following lists:
    [sign, int deg, frac deg]
    [sign, int deg, int min, frac min]
    [sign, int deg, int min, int sec, frac sec]
    where:
        all values are strings
        sign is one of ('', '+' or '-')
        frac <whatever> includes a leading decimal point
    
    error conditions:
        raises ValueError if the string cannot be parsed
    """
    assert isinstance(dmsStr, str)
    m = _DegRE.match(dmsStr) or _DegMinRE.match(dmsStr) or _DegMinSecRE.match(dmsStr)
    if m == None:
        raise ValueError, "splitDMSStr cannot parse %s as a sexagesimal string" % (dmsStr)
    matchSet = list(m.groups())
    if matchSet[-1] == None:
        matchSet[-1] = ''
    return matchSet

_FloatRE = re.compile(r'^\s*[-+]?[0-9]*\.?[0-9]*(e[-+]?)?[0-9]*\s*$', re.IGNORECASE)
_FloatNoExpRE = re.compile(r'^\s*[-+]?[0-9]*\.?[0-9]*\s*$')
def floatFromStr(astr, allowExp=1):
    """Convert a string representation of a number to a float;
    unlike float(), partial representations (such as "", "-", "-.e") are taken as 0
    and "nan" is forbidden.

    error conditions:
        raises ValueError if astr cannot be converted
    """
    if allowExp:
        match = _FloatRE.match(astr)
    else:
        match = _FloatNoExpRE.match(astr)
    
    
    if match == None:
        raise ValueError, "cannot convert :%s: to a float" % (astr)
        
    try:
        return float(astr)
    except:
        # partial float
        return 0.0

_IntRE = re.compile(r'^\s*[-+]?[0-9]*\s*$')
def intFromStr(astr):
    """Convert a string representation of a number to an integer;
    unlike int(), the blank string and "+" and "-" are treated as 0

    error conditions:
        raises ValueError if astr cannot be converted
    """
    if _IntRE.match(astr) == None:
        raise ValueError, "cannot convert :%s: to an integer" % (astr)

    try:
        return int(astr)
    except:
        # partial int
        return 0

def quoteStr(astr, escChar='\\', quoteChar='"'):
    """Escape all instances of quoteChar and escChar in astr
    with a preceding escChar and surrounds the result with quoteChar.
    
    Examples:
    astr = 'foo" \bar'
    quoteStr(astr) = '"foo\" \\bar"'
    quoteStr(astr, escChar = '"') = '"foo"" \bar"'

    This prepares a string for output.
    """
    if escChar != quoteChar:
        # escape escChar
        astr = astr.replace(escChar, escChar + escChar)
    # escape quoteChar and surround the result in quoteChar
    return quoteChar + astr.replace(quoteChar, escChar + quoteChar) + quoteChar

def unquoteStr(astr, escChar='\\', quoteChars='"\''):
    """Remove quotes from a string and unescapes contained escaped quotes.
    
    Based on email.unquote.
    """
    if len(astr) > 1:
        for quoteChar in quoteChars:
            if astr.startswith(quoteChar) and astr.endswith(quoteChar):
                return astr[1:-1].replace(escChar + escChar, escChar).replace(escChar + quoteChar, quoteChar)
    return astr

def strFromException(exc):
    """Unicode-safe replacement for str(exception)"""
    try:
        return str(exc)
    except Exception:
        try:
            return ",".join([unicode(s) for s in exc.args])
        except Exception:
            # in case exc is some unexpected type
            return repr(exc)
      
def _getDMSFields (decDeg, nFields=3, precision=1):
    """Return a string representation of dms fields for decDeg.

    Inputs:
    - decDeg: value in decimal degrees or hours
    - nFields: number of fields; must be >= 1 (and >3 is probably never used)
    - precision: number of digits after the decimal point in the last field;
        if 0, no decimal point is printed; must be >= 0
     
    Returns:
    - signStr: "" or "-"
    - fieldStrs: string value of each field (all positive)
    
    To compute dms a a string, use: signStr + ":".join(fieldStrs)
    This routine doesn't take that step to allow omitting fields
    whose value is zero (e.g. see dmsStrFromDeg and dmsStrFromSec).
    
    Error conditions:
    - Raises ValueError if precision < 0
    """
    if nFields < 1:
        raise ValueError("nFields=%r; must be >= 1" % (nFields,))
    if precision < 0:
        raise ValueError("precision=%r; must be >= 0" % (precision,))

    if decDeg < 0:
        signStr = "-"
    else:
        signStr = ""

    if nFields == 1:
        retNum = round(abs(decDeg), precision)
        retStr = "%.*f" % (precision, retNum)
        return signStr, [retStr]
    
    # compute a list of output fields; all but the last one are integer
    remVal = abs(decDeg)
    fieldNums = []
    for fieldNum in range(nFields-1):
        (intVal, remVal) = divmod (abs(remVal), 1.0)
        intVal = int(intVal)
        fieldNums.append(intVal)
        remVal *= 60.0
    fieldNums.append(round(remVal, precision))
    
    # handle overflow
    incrPrevField = False
    for ind in range(nFields-1, -1, -1):
        if incrPrevField:
            fieldNums[ind] += 1
        if (ind > 0) and (fieldNums[ind] >= 60):
            fieldNums[ind] -= 60
            incrPrevField = True
        else:
            incrPrevField = False

    # compute fieldStrs
    if precision > 0:
        minFloatWidth = precision + 3
    else:
        minFloatWidth = 2
    fieldStrs = ["%d" % (fieldNums[0],)]
    for numVal in fieldNums[1:-1]:
        fieldStrs.append("%02d" % (numVal,))
    fieldStrs.append("%0*.*f" % (minFloatWidth, precision, fieldNums[-1]))
            
    return signStr, fieldStrs

def _assertTest():
    """Run a test by comparing results to those expected and only failing if something is wrong.
    
    Use _printTest to generate data for this test.
    """
    # format is a set of lists of:
    # - dms string
    # - comment
    # - should work (True or False)
    # - splitStr: the expected output from splitDMSStr
    # - degVal: the expected output from degFromDMSStr
    # - secVal: the expected output from secFromDMSStr
    # - neatStr: the expected output from neatenDMSStr
    # - a list of three expected outputs from dmsStrFromDeg with nFields=3 and:
    #   - precision = 0
    #   - precision = 1
    #   - precision = 2
    testSet = (
        ['::', '', True, [''    , '', '', '', ''], 0.0, 0.0, '0:00:00', ['0:00:00', '0:00:00.0', '0:00:00.00']],
        ['-::', '', True, ['-', '', '', '', ''], -0.0, -0.0, '0:00:00', ['0:00:00', '0:00:00.0', '0:00:00.00']],
        ['-0:00:00.01', '', True, ['-', '0', '00', '00', '.01'], -2.7777777777777775e-06, -0.01, '-0:00:00.01', ['-0:00:00', '-0:00:00.0', '-0:00:00.01']],
        [' +1', '', True, ['+', '1', ''], 1.0, 1.0, '1', ['1:00:00', '1:00:00.0', '1:00:00.00']],
        ['-1.2345', '', True, ['-', '1', '.2345'], -1.2344999999999999, -1.2344999999999999, '-1.2345', ['-1:14:04', '-1:14:04.2', '-1:14:04.20']],
        ['-123::', '', True, ['-', '123', '', '', ''], -123.0, -442800.0, '-123:00:00', ['-123:00:00', '-123:00:00.0', '-123:00:00.00']],
        ['-123:4', 'make sure seconds field is not 60 from dmsStrFromDeg', True, ['-', '123', '4', ''], -123.06666666666666, -7384.0, '-123:04', ['-123:04:00', '-123:04:00.0', '-123:04:00.00']],
        ['-123:45', '', True, ['-', '123', '45', ''], -123.75, -7425.0, '-123:45', ['-123:45:00', '-123:45:00.0', '-123:45:00.00']],
        ['-123:4.56789', '', True, ['-', '123', '4', '.56789'], -123.0761315, -7384.5678900000003, '-123:04.56789', ['-123:04:34', '-123:04:34.1', '-123:04:34.07']],
        ['-123:45.6789', '', True, ['-', '123', '45', '.6789'], -123.761315, -7425.6788999999999, '-123:45.6789', ['-123:45:41', '-123:45:40.7', '-123:45:40.73']],
        ['1:2:', '', True, ['', '1', '2', '', ''], 1.0333333333333334, 3720.0, '1:02:00', ['1:02:00', '1:02:00.0', '1:02:00.00']],
        ['1:2:3', '', True, ['', '1', '2', '3', ''], 1.0341666666666667, 3723.0, '1:02:03', ['1:02:03', '1:02:03.0', '1:02:03.00']],
        ['1:2:3.456789', '', True, ['', '1', '2', '3', '.456789'], 1.0342935525000001, 3723.4567889999998, '1:02:03.456789', ['1:02:03', '1:02:03.5', '1:02:03.46']],
        ['1:23:4', '', True, ['', '1', '23', '4', ''], 1.3844444444444444, 4984.0, '1:23:04', ['1:23:04', '1:23:04.0', '1:23:04.00']],
        ['1:23:45', '', True, ['', '1', '23', '45', ''], 1.3958333333333333, 5025.0, '1:23:45', ['1:23:45', '1:23:45.0', '1:23:45.00']],
        ['123:45:6.789', '', True, ['', '123', '45', '6', '.789'], 123.75188583333333, 445506.78899999999, '123:45:06.789', ['123:45:07', '123:45:06.8', '123:45:06.79']],
        ['123:45:56.789', '', True, ['', '123', '45', '56', '.789'], 123.76577472222222, 445556.78899999999, '123:45:56.789', ['123:45:57', '123:45:56.8', '123:45:56.79']],
        ['-0::12.34', 'bug test; the sign must be retained', True, ['-', '0', '', '12', '.34'], -0.0034277777777777779, -12.34, '-0:00:12.34', ['-0:00:12', '-0:00:12.3', '-0:00:12.34']],
        ['-::12.34', 'a weird gray area, but it works', True, ['-', '', '', '12', '.34'], -0.0034277777777777779, -12.34, '-0:00:12.34', ['-0:00:12', '-0:00:12.3', '-0:00:12.34']],
        ['::12.34', '', True, ['', '', '', '12', '.34'], 0.0034277777777777779, 12.34, '0:00:12.34', ['0:00:12', '0:00:12.3', '0:00:12.34']],
        ['1:23.4567', '', True, ['', '1', '23', '.4567'], 1.3909450000000001, 83.456699999999998, '1:23.4567', ['1:23:27', '1:23:27.4', '1:23:27.40']],
        ['-1.234567', '', True, ['-', '1', '.234567'], -1.234567, -1.234567, '-1.234567', ['-1:14:04', '-1:14:04.4', '-1:14:04.44']],
        ['-1:abadstr', 'invalid characters', False, None, None, None, None, None],
        ['-1:2343:24', 'too many minutes digits', False, None, None, None, None, None],
        ['1:-1:24', 'minus sign in wrong place', False, None, None, None, None, None],
    )
    def locAssert(fmt, res, func, *args, **kargs):
        assert fmt % (res,) == fmt % (func(*args, **kargs),), "%r != %r = %s(*%r, **%r)" % (res, func(*args, **kargs), func.__name__, args, kargs)
    
    nErrors = 0
    for testStr, commentStr, isOK, splitStr, degVal, secVal, neatStr, dmsStr02 in testSet:
        try:
            locAssert("%r", splitStr, splitDMSStr, testStr)
            locAssert("%.8g", degVal, degFromDMSStr, testStr)
            locAssert("%.8g", secVal, secFromDMSStr, testStr)
            locAssert("%r", neatStr, neatenDMSStr, testStr)
            locAssert("%r", dmsStr02[0], dmsStrFromDeg, degVal, 3, 0)
            locAssert("%r", dmsStr02[1], dmsStrFromDeg, degVal, 3, 1)
            locAssert("%r", dmsStr02[2], dmsStrFromDeg, degVal, 3, 2)
            if not isOK:
                print "unexpected success on %r" % testStr
                nErrors += 1
        except StandardError, e:
            if isOK:
                raise
                print "unexpected failure on %r\n\t%s\nskipping other tests on this value" % (testStr, e)
                nErrors += 1
    
    if nErrors == 0:
        print "RO.StringUtil passed"
    else:
        print "RO.StringUtil failed with %d errors" % nErrors
            

def _printTest(dmsSet = None):
    """Print the results of running each routine on a set of test data.
    Data format is a list of tuples, each containing two elements:
        dmsStr to test, a comment
        
    The output is in the format used by _assertTest, but please use this with great caution.
    You must examine the output very carefully to confirm it is correct before updating _assertTest!
    """
    print "Exercising RO string utilities"
    if not dmsSet:
        dmsSet = (
            ("::", ""),
            ("-::", ""),
            ('-0:00:00.01', ""),
            (" +1", ""),
            ('-1.2345', ''),
            ('-123::', ''),
            ('-123:4', 'make sure seconds field is not 60 from dmsStrFromDeg'),
            ('-123:45', ''),
            ('-123:4.56789', ''),
            ('-123:45.6789', ''),
            ('1:2:', ''),
            ('1:2:3', ''),
            ('1:2:3.456789', ''),
            ('1:23:4', ''),
            ('1:23:45', ''),
            ('123:45:6.789', ''),
            ('123:45:56.789', ''),
            ('-0::12.34', 'bug test; the sign must be retained'),
            ('-::12.34', 'a weird gray area, but it works'),
            ('::12.34', ''),
            ('1:23.4567', ''),
            ('-1.234567', ''),
            ('-1:abadstr', 'invalid characters'),
            ('-1:2343:24', 'too many minutes digits'),
            ('1:-1:24', 'minus sign in wrong place'),
        )

    for testStr, commentStr in dmsSet:
        # note: if splitDMSStr succeeds, then the other calls all should succeed
        if checkDMSStr(testStr):
            try:
                itemList = splitDMSStr(testStr)
                deg = degFromDMSStr (testStr)
                sec = secFromDMSStr (testStr)
                neatStr = neatenDMSStr (testStr)
                outDMSStr = []
                for prec in range(3):
                    outDMSStr.append(dmsStrFromDeg(deg, precision=prec))
                print "[%r, %r, True, %r, %r, %r, %r, %r]," % (testStr, commentStr, itemList, deg, sec, neatStr, outDMSStr)
            except StandardError, e:
                print "unexpected failure on %r (%s); error = %s" % (testStr, commentStr, e)
        else:
            print "[%r, %r, False, %r, %r, %r, %r, %r]," % tuple([testStr, commentStr] + [None]*5)

if __name__ == "__main__":
    doPrint = False
    if doPrint:
        _printTest()
    else:
        _assertTest()
