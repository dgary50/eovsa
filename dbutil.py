'''
   Module for creating and manipulating SQL database queries'''
#
# History:
#   2014-Nov-29  DG
#      New module.
#   2015-Feb-20  DG
#      Changed hard-coded v38 to v39 in table name.  This requirement has to be
#      removed by new code on the SQL server, planned but not yet implemented.
#   2015-Apr-01  DG
#      At some point, get_dbrecs acquired a version keyword to avoid hard-coded
#      version number.  Added do_query() function for general query.
#   2015-Apr-02  DG
#      Finally figured out a way to get the right table version for a given
#      timestamp.  Added routine find_table_version() to accomplish it.
#   2016-Aug-03  DG
#      Slight change in get_dbrecs(), to call find_table_version(), if the
#      version is not given.
#   2016-Aug-04  DG
#      Another change to allow get_dbrecs() to be called with a Time()
#      object or even a timerange (which means nrecs need not be given).
#   2016-Aug-06  DG
#      Made get_chi() and a14_wsram() version independent.
#   2016-Nov-26  DG
#      Added get_motor_currents().
#   2017-Apr-27  DG
#      Added return of average wind speed to a14_wscram()
#   2017-May-16  DG
#      Added get_reboot() for finding ROACH reboot times
#   2017-Aug-06  DG
#      Changed get_dbrecs() so that timerange is inclusive, i.e.
#      returns data for both start and end second.
#

import stateframedef
import util
from util import Time

def get_cursor():
    ''' Connect to the SQL database and return a cursor for access to it
    '''
    cnxn = stateframedef.pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;")
    return cnxn.cursor()
    
def find_table_version(cursor,timestamp,scan_header=False):
    ''' Searches dimension-1 tables for all versions in the database
        to find the one containing the given timestamp.  Returns the
        version number as a string, e.g. '51'
    '''
    import fnmatch
    filtstr = 'fV??_vD1'
    if scan_header:
        filtstr = 'hV??_vD1'
    cursor.tables(tableType='VIEW')
    rows = cursor.fetchall()
    version = None
    for row in rows:
        tbl = row[2]
        if fnmatch.filter([tbl],filtstr) != []:
            # This is a "version" dimension-1 table, so get its start time
            try:
                cursor.execute('select top 1 Timestamp from '+tbl)
                tstamp = cursor.fetchone()[0]
                if tstamp < timestamp:
                    mytbl = tbl
                    version = mytbl[2:4]
            except:
                pass
    return version
    
def get_dbrecs(cursor=None,version=None,dimension=None,timestamp=None,nrecs=None):
    ''' Fairly general routine for fetching a contiguous block of data and returning
        it as a dictionary of arrays of size nrecs x dimension.
        
        Note: timestamp can be given as a single LabVIEW timestamp, or a
        single Time() object, or as a two-element Time() object representing
        a timerange.  If the latter, nrecs is determined from the timerange.
    '''
    te = None
    if type(timestamp) == util.Time:
        try:
            if len(timestamp) == 2:
                # This is a timerange as Time object.  Generate nrecs from time difference (in s)
                ts = timestamp[0].lv
                nrecs = int(round(timestamp[1].lv - timestamp[0].lv)) + 1
            else:
                print 'Too many times in Time() object.'
                return {}
        except:
            # This is a single Time object
            ts = timestamp.lv
    else:
        ts = timestamp
    if type(cursor) != stateframedef.pyodbc.Cursor:
        print 'No database open'
        return {}
    if ts is None:
        print 'A timestamp must be given.'
        return {}
    if version is None:
        version = int(find_table_version(cursor,ts))
    if type(version) != int:
        print 'Version must be int type.'
        return {}
    if type(dimension) != int:
        print 'Dimension must be int type.'
        return {}
    if type(nrecs) != int:
        print 'NRecs must be int type.'
        return {}
    nvals = dimension*nrecs
    # Generate table name
    table = 'fV'+str(version)+'_vD'+str(dimension)
    # Generate query
    query = 'select top '+str(nvals)+' * from '+table+' where timestamp >= '+str(ts)
    try:
        cursor.execute(query)
    except:
        print 'Query',query.upper(),'returned an error.'
        print stateframedef.sys.exc_info()[0]
    # Extract the data
    data = stateframedef.numpy.transpose(stateframedef.numpy.array(cursor.fetchall(),'object'))
    # Override nrecs with the number of records actually read (could be less than requested)
    nrecs = len(data[0])/dimension
    # Get names from description
    names = stateframedef.numpy.array(cursor.description)[:,0]
    # Reshape data array for zipping into dictionary.  Each dictionary entry will be
    # an array of size nrecs x dimension.
    if dimension > 1:
        data.shape = (len(names),nrecs,dimension)
    else:
        data.shape = (len(names),nrecs)
    # Create the dictionary
    outdict = dict(zip(names,data))
    return outdict
    
def do_query(cursor,query):
    ''' Executes the supplied query on an already open database pointed
        to by cursor.  Returns the result of the query as a dictionary
        (could be an empty dictionary if no results were returned).
        Also returns a message indicating success or an error:
        
         outdict, msg = do_query(cursor, query) 
    '''
    import sys
    try:
        cursor.execute(query)
        data = stateframedef.numpy.transpose(stateframedef.numpy.array(cursor.fetchall(),dtype='object'))
        names = stateframedef.numpy.array(cursor.description)[:,0]
        result = dict(zip(names,data))
        msg = 'Success'
    except:
        result = {}
        msg = 'Error: '+str(sys.exc_info()[1])
    return result,msg
    
def a14_wscram(trange):
    ''' Get the Antenna 14 windscram state, and the average wind speed, for a 
        given time range.
        
        Returns:
           times      as a Time() object, or error message if failure
           wscram     array of windscram state, 0 = not in wind scram, 1 = in windscram
           avgwind    array of average wind speeds, in MPH, or error message if failure
    '''
    tstart,tend = [str(i) for i in trange.lv]
    cursor = get_cursor()
    ver = find_table_version(cursor,trange[0].lv)
    query = 'select Timestamp,Ante_Fron_Wind_State from fV'+ver+'_vD15 where (I15 = 13) and Timestamp between '+tstart+' and '+tend
    data, msg = do_query(cursor, query)
    if msg == 'Success':
        times = Time(data['Timestamp'].astype('int'),format='lv')
        wscram = data['Ante_Fron_Wind_State']
    else:
        return 'Error: '+msg, None, None
    query = 'select Timestamp,Sche_Data_Weat_AvgWind from fV'+ver+'_vD1 where Timestamp between '+tstart+' and '+tend
    data, msg = do_query(cursor, query)
    if msg == 'Success':
        avgwind = data['Sche_Data_Weat_AvgWind']
    else:
        return times,wscram,'Error: '+msg
    cursor.close()
    return times,wscram,avgwind
    
def get_chi(trange):
    ''' Get the parallactic angle for all antennas (ntimes x 16) for a
        given time range (returns times and parallactice angle--radians)
    '''
    tstart,tend = [str(i) for i in trange.lv]
    cursor = get_cursor()
    ver = find_table_version(cursor,trange[0].lv)
    query = 'select Timestamp,I16,Sche_Data_Chi from fV'+ver+'_vD16 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
    data, msg = do_query(cursor, query)
    cursor.close()
    if msg == 'Success':
        times = Time(data['Timestamp'].astype('int'),format='lv')[::16]
        chi = data['Sche_Data_Chi']
        nt = len(chi)/16
        chi.shape = (nt,16)
        return times,chi

def get_motor_currents(trange):
    ''' Get the Azimuth and Elevation motor currents for all antennas (ntimes x 15) for a
        given time range (returns times, azimuth motor current, and elevation motor current)
    '''
    tstart,tend = [str(i) for i in trange.lv]
    cursor = get_cursor()
    ver = find_table_version(cursor,trange[0].lv)
    query = 'select Timestamp,Ante_Cont_AzimuthMotorCurrent,Ante_Cont_ElevationMotorCurrent from fV'+ver+'_vD15 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
    data, msg = do_query(cursor, query)
    cursor.close()
    if msg == 'Success':
        times = Time(data['Timestamp'].astype('int'),format='lv')[::15]
        az = data['Ante_Cont_AzimuthMotorCurrent']
        el = data['Ante_Cont_ElevationMotorCurren']
        nt = len(az)/15
        az.shape = (nt,15)
        el.shape = (nt,15)
        return times,az,el
    else:
        print msg
        return None,None,None

def get_reboot(trange,previous=False):
    ''' Get the times of any ROACH (correlator) reboots in the given timerange
    
        Returns Time() object of reboots, or None. If previous is True, it
        returns the time of the previous reboot.
    '''
    import numpy as np
    t0, t1 = trange.lv.astype(np.int)
    tmjd = trange[0].mjd
    cursor = get_cursor()
    ver = find_table_version(cursor,t0,scan_header=True)
    query = 'select Timestamp,TimeAtAcc0 from hV'+ver+'_vD1 where Timestamp between '+str(t0)+' and '+str(t1)+' order by Timestamp'
    data, msg = do_query(cursor, query)
    t0, idx = np.unique(data['TimeAtAcc0'],return_index=True)
    t_reboot = data['TimeAtAcc0'][idx].astype(float)
    if previous:
        pass
    else:
        if t_reboot.size == 1:
            if t_reboot < tmjd:
                return None
        elif t_reboot[0] < tmjd:
            t_reboot = t_reboot[1:]
    return Time(t_reboot,format='mjd')

