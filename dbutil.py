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

import stateframedef

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
    '''
    if type(cursor) != stateframedef.pyodbc.Cursor:
        print 'No database open'
        return {}
    if type(version) != int:
        print 'Version must be int type.'
        return {}
    if type(dimension) != int:
        print 'Dimension must be int type.'
        return {}
    if type(nrecs) != int:
        print 'NRecs must be int type.'
        return {}
    if type(timestamp) is None:
        print 'A timestamp must be given.'
        return {}
    nvals = dimension*nrecs
    # Generate table name
    table = 'fV'+str(version)+'_vD'+str(dimension)
    # Generate query
    query = 'select top '+str(nvals)+' * from '+table+' where timestamp >= '+str(timestamp)
    try:
        cursor.execute(query)
    except:
        print 'Query',query.upper(),'returned an error.'
        print stateframedef.sys.exc_info()[0]
    # Extract the data
    data = stateframedef.numpy.transpose(stateframedef.numpy.array(cursor.fetchall(),'object'))
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
    
