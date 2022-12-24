#
#  Routines for defining stateframe and scan header SQL database representations
#
#  2015-May-29  DG
#    Converted from using datime() to using Time() based on astropy.
#  2015-Oct-17  DG
#    Latest stateframe (v60) revealed a bug in logic of startbyte().  Also
#    fixed minor glitch is log2sql(), which was looking for stateframe.dt,
#    which no longer exists.  Instead, this routine now imports datetime 
#    directly.
#  2016-Apr-02  DG
#    The v64 stateframe caused problems due to putting a string (byte array)
#    inside an array.  The problem became really hard to fix due to an
#    unrelated bug in the code in sfdef().  After fixing that bug, a
#    special case for Parser Command string had to be added that actually
#    overrides the stateframe dictionary to include the string size as a
#    separate item.  Now I have to go back and recreate the v64 table, and
#    upload all of the stateframe data recorded in the log files for the
#    past month...
#  2016-Apr-05  DG
#    Added new routine drop_deftable(), to drop an existing table definition
#    and delete all tables associated with it.  This is a dangerous command,
#    so it requests confirmation from the user via the keyboard.  Also fixed
#    bug in sfdef() that occurred when the passed-in dictionary is a scan header.
#

import stateframe, struct, numpy, pyodbc, datetime
import read_xml2 as rxml
import sys, os, time, glob

#=============== get_start_byte ===============
def get_start_byte(c):
    '''Find the start byte of the first endpoint in dict c
    '''
    if not (type(c) is dict):
        print 'Argument must be a dictionary'
        return -1
    keys = c.keys()
    for key in keys:
        if type(c[key]) is dict:
            sbyte = get_start_byte(c[key])
            break
        else:
            sbyte = c[key][1]
            break
    return sbyte

#=============== walk keys ===============
def walk_keys(c,inkey,indim=1,pre='',nbytes=0):
    ''' Handle an entity c
        If a 'dict', recursively handle its keys
        If a 'list', check if it is an endpoint or an array
    '''
    # Remove any white spaces in inkey
    inkey = inkey.replace(' ','')
    # Remove any '.' in inkey but capitalize word after '.'
    if inkey.find('.') != -1:
        # Although the statement below works even when no dots, it can
        # lower the case of some characters in the string, which is not wanted
        inkey = inkey.replace('.',' ').title().replace(' ','')
    dtype_dict = {'s':['string',1],'B':['byte',1],'H':['u16',2],'I':['u32',4],
                  'h':['i16',2],'i':['i32',4],'f':['float',4],'d':['double',8]}
    dict_keys = ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']

    mylist = []
    if type(c) is dict:
        keys = c.keys()
        for key in keys:
            #print '+'+inkey+'_'+key
            newlist = walk_keys(c[key],key,indim,pre+key[0:4]+'_',nbytes)
            if newlist is None:
                pass
            else:
                mylist += newlist
    else:
        # This is a list, so either it is an endpoint or an array
        if type(c[0]) is dict:
            # This is an array of dicts--no extension of pre string
            indim = len(c)
            nbytes = get_start_byte(c[1]) - get_start_byte(c[0])
            newlist = walk_keys(c[0],inkey,indim,pre,nbytes)
            if newlist is None:
                pass
            else:
                mylist += newlist
        else:
            # This must be an endpoint
            # Check if the endpoint is itself an array
            # Prefix contains the current key, so pare it off -- p is the actual prefix
            p = pre[:-min(5,len(inkey)+1)]
            if len(c) == 3:
                # This is an array of values.  
                if len(c[2]) == 1:
                    # This is a 1D array
                    indim = c[2][0]
                    if c[0][-1] == 's':
                        nbyt = int(c[0][:-1])
                        dtype = 'string'
                    else:
                        dtype,nbyt = dtype_dict[c[0][-1]]
                    # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                    if nbytes == 0:
                        mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1],nbyt,c[0],p+inkey]))]
                    else:
                        mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1],nbytes,c[0],p+inkey]))]
                else:
                    # This is a 2D array. Unfortunately, it is here
                    # that there are several special-case exceptions to
                    # deal with. These have to be handled explicitly.
                    if inkey == 'UVW':
                        d1 = c[2][0]  # first dimension
                        d2 = c[2][1]  # second dimension
                        indim = d2
                        nbytes = struct.calcsize(str(d2)+c[0][-1])
                        # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                        for i in range(d1):
                            dtype,nbyt = dtype_dict[c[0][-1]]
                            mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1]+nbytes*i,nbyt,c[0],p+inkey[i]]))]
                    elif inkey == 'Delay':
                        d2 = c[2][0]  # first dimension
                        d1 = c[2][1]  # second dimension
                        indim = d2
                        nbytes = struct.calcsize(str(d2)+c[0][-1])
                        # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                        for i in range(d1):
                            dtype,nbyt = dtype_dict[c[0][-1]]
                            mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1]+nbyt*i,nbyt*2,c[0][-1],p+inkey+str(i)]))]
                    elif inkey == 'Antpos':
                        d1 = c[2][0]  # first dimension
                        d2 = c[2][1]  # second dimension
                        indim = d2
                        nbytes = struct.calcsize(str(d2)+c[0][-1])
                        ants = ['X','Y','Z']
                        # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                        for i in range(d1):
                            dtype,nbyt = dtype_dict[c[0][-1]]
                            mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1]+nbyt*i,nbyt*d1,c[0][-1],p+inkey+ants[i]]))]
                    elif inkey == 'Ephem' or inkey == 'SatEphem':
                        d1 = c[2][0]  # first dimension
                        d2 = c[2][1]  # second dimension
                        indim = d2
                        nbytes = struct.calcsize(str(d2)+c[0][-1])
                        num = ['Time','RA','Dec']
                        # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                        for i in range(d1):
                            dtype,nbyt = dtype_dict[c[0][-1]]
                            mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1]+nbyt*i,nbyt*d1,c[0][-1],p+inkey+num[i]]))]
            else:
                # This is a single endpoint
                if c[0][-1] == 's':
                    nbyt = int(c[0][:-1])
                    dtype = 'string'
                else:
                    dtype,nbyt = dtype_dict[c[0][-1]]
                # ['dimension','datatype','fieldbytes','startbyte','dimoffset','pytype','fieldname']
                if nbytes == 0:
                    mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1],nbyt,c[0],p+inkey]))]
                else:
                    mylist += [dict(zip(dict_keys,[indim,dtype,nbyt,c[1],nbytes,c[0],p+inkey]))]
    return mylist
    
#=============== sfdef ===============
def sfdef(sf=None):
    ''' Main routine to generate a dictionary list representing the
        stateframedef table (outlist) for either a stateframe or a 
        scan_header dictionary, and a dictionary list of rearrangements 
        to the binary data, brange.  If sf is omitted, the current 
        stateframe.xml file is used generate it.
    '''
    # helper routines list_sort, sub_list, and rd_sfdef
    def list_sort(mylist,key):
        ''' Given a key whose value is numerical, sort the list
            of dictionaries mylist according to the key and return
            the new list of dictionaries
        '''
        val = []
        newlist = []
        for d in mylist:
            val.append(d[key])
        idx = numpy.array(val).argsort()
        for i in idx:
            newlist += [mylist[i]]
        return newlist

    def sub_list(mylist):
        ''' Given a sorted list of dictionaries mylist, separate them
            into lists of a given dimension.  NB: entries with dimension
            greater than 50 are ignored and will NOT appear in the RDBMS
        '''
        dim = mylist[0]['dimension']
        lists = []
        sublist = []
        for d in mylist:
            if d['dimension'] == dim:
                sublist += [d]
            elif d['dimension'] > 50:
                # Skip any entries with dimension > 50
                pass
            else:
                lists += [sublist]
                sublist = [d]
                dim = d['dimension']
        lists += [sublist]            
        return lists

    def rd_sfdef(sf):
        ''' Create initial list of dictionaries representing each variable
            in the given dictionary sf.  This routine works with either
            a stateframe dictionary or a scan_header dictionary.  However,
            if sf is not given, sf is obtained from the current stateframe
            version, and that only works for the stateframe.  
        '''
        mylist = walk_keys(sf,'Stateframe')
        return mylist

    if sf is None:
        accini = stateframe.rd_ACCfile()
        sf = accini['sf']

    # This is something needed if the dictionary is a stateframe--will fail
    # gracefully and skip this if the dictionary is a scan header
    try:
        # This is a major breaking of the scheme, but the Parser string in the Antenna
        # array requires special handling:
        for ant in range(15):
            # For each antenna, add an entry to the stateframe representing the Command string size
            # located 4 bytes before the Command string
            t,loc = sf['Antenna'][ant]['Parser']['Command']                     # Get Command string location
            sf['Antenna'][ant]['Parser'].update({'Command_Size':['i',loc-4]})   # Insert another entry pointing to 4 bytes before
    except KeyError:
        # The dictionary is probably a scan header
        pass
    # First parse the XML file to get a list of dictionaries for
    # stateframedef SQL table
    mylist = rd_sfdef(sf)

    # Sort the list of dictionaries by dimension
    newlist = list_sort(mylist,'dimension')
    # Separate sorted list into separate lists according to dimension 
    lists = sub_list(newlist)
    brange = []
    outlist = []
    # Loop over separated lists
    for alist in lists:
        # Sort the list according to startbyte
        slist = list_sort(alist,'startbyte')
        outlist += slist
        # Step through ordered list to look for gaps, and make map of
        # start and end bytes for contiguous byte ranges
        sbyte = slist[0]['startbyte']
        for i,item in enumerate(slist[:-1]):
            if len(item['pytype']) > 1 and item['datatype'] != 'string':
                # Case of "in-place" non-string array
                ebyte = item['startbyte']+item['fieldbytes']*item['dimension']
            else:
                ebyte = item['startbyte']+item['fieldbytes']
            if ebyte != slist[i+1]['startbyte']:
                # There is a gap between this and next item, but this might
                # be due to dimensionality of this section.  Add bytes for
                # dimensionality and check again.
                if len(item['pytype']) > 1 and item['datatype'] != 'string':
                    # Case of "in-place" non-string array -- do nothing
                    pass
                else:
                    ebyte = (ebyte - sbyte)*item['dimension'] + sbyte
                if ebyte != slist[i+1]['startbyte']:
                    # This is indeed a gap, so save byte range and start a
                    # new byte range
                    brange.append({'sbyte':sbyte,'ebyte':ebyte})
                    sbyte = slist[i+1]['startbyte']
        item = slist[-1]
        ebyte = item['startbyte']+item['fieldbytes']
        # Take care of dimensionality of last item
        if len(item['pytype']) > 1:
            ebyte = (ebyte - item['startbyte'])*item['dimension'] + item['startbyte']
            #print '5:',sbyte, (ebyte-slist[-1]['startbyte'])/slist[-1]['dimension'] + slist[-1]['startbyte'], slist[-1]['dimension']
        else:
            ebyte = (ebyte - sbyte)*item['dimension'] + sbyte
            #print '5:',sbyte, (ebyte-sbyte)/slist[-1]['dimension'] + sbyte, slist[-1]['dimension']
        brange.append({'sbyte':sbyte,'ebyte':ebyte})
    # At this point, we should have a complete map of the ordered byte
    # ranges into original binary data, in variable brange.
    return brange, outlist

#=============== startbyte ===============
def startbyte(outlist):
    # The input startbyte has to be transferred to origbyte, and
    # we have to do another pass through to calculate output startbyte
    curpos = 0
    dim = outlist[0]['dimension']
    offset = outlist[0]['dimoffset']
    outlist[0].update({'origbyte':outlist[0]['startbyte']})
    outlist[0].update({'startbyte':curpos})
    curpos += outlist[0]['fieldbytes']
    for i in range(1,len(outlist)):
        outlist[i].update({'origbyte':outlist[i]['startbyte']})
        outlist[i].update({'startbyte':curpos})
        if outlist[i]['dimension'] == 1:
            # For dimension 1 items, dimoffset is irrelevant, so set to 0
            outlist[i]['dimoffset'] = 0
        # Conditions that can mark the end of a block of data
        cond1 = outlist[i]['dimension'] != dim    # current line has a dimension change
        cond2 = outlist[i]['dimoffset'] != offset # current line has a dimoffset change
        if outlist[i]['datatype'] == 'string':
            # Requires special handling, because this is formally an array, with an array
            # dimension in the data, but it does not signify the end of a block of data
            cond3 = ((outlist[i]['origbyte']  - outlist[i-1]['origbyte']) > 
                     (outlist[i]['startbyte'] - outlist[i-1]['startbyte'] + 4))  # Gap in origbyte (accounting for string dimension)
            cond4 = len(outlist[i]['pytype']) > 1 and outlist[i]['datatype'] != 'string'     # A non-string array
        else:
            cond3 = ((outlist[i]['origbyte']  - outlist[i-1]['origbyte']) > 
                     (outlist[i]['startbyte'] - outlist[i-1]['startbyte']))  # Gap in origbyte
            cond4 = len(outlist[i]['pytype']) > 1 # An array
        if (cond1 or cond2 or cond3 or cond4):
            # Update the current position by n-1 x offset of previous line
            curpos += outlist[i-1]['dimoffset']*(dim-1)
            outlist[i].update({'startbyte':curpos})
            dim = outlist[i]['dimension']
            offset = outlist[i]['dimoffset']
        curpos += outlist[i]['fieldbytes']

#=============== outlist2table ===============
def outlist2table(outlist,version):
    ''' Convert the completed outlist to an SQL table, which is a list of the stateframedef 
        or scanheaderdef table parameters, used by load_deftable().  
    '''
    idx = 1
    dim = 1
    # If version is an integer (i.e. no float remainder) then write it out as
    # an integer (this converts, e.g., 25.0 to 25, but leaves 25.1 as is)
    ver = version
    tbl = []
    if (version % 1) == 0:
        ver = int(version)
    for i,item in enumerate(outlist):
        if item['dimension'] != dim:
            # This generates the "fieldnum," which is just the ordinal position in a table
            # of a given dimension.  It resets to 1 whenever a dimension change occurs in
            # the ordered table.
            idx = 1
            dim = item['dimension']
        # print '{:8} {:8} {:8} {:8} {:8} {:3} {:8} {:8} {:8}    {:}'.format(i+1,0,ver,item['dimension'],item['datatype'],item['fieldbytes'],item['dimoffset'],item['startbyte']+1,idx,item['fieldname'])
        tbl.append([ver,item['dimension'],item['datatype'],item['fieldbytes'],item['dimoffset'],item['startbyte']+1,idx,item['fieldname']])
        idx += 1
    return tbl
    
#=============== transmogrify ===============
def transmogrify(indata,brange):
    ''' Rearrange incoming stateframe binary data to order it
        in the way described by the stateframedef table.
           indata: the incoming stateframe binary data
           brange: the list of ranges, as a list of dictionaries

        Returns outdata, the rearranged binary data
    '''
    outdata = ''
    for item in brange:
        outdata += indata[item['sbyte']:item['ebyte']]

    return outdata

#=============== old_version_test ===============
def old_version_test(sflog=None,sfxml=None,outbinfile=None,outtabfile=None):
    ''' Read stateframe log files of older versions and
        create output file of rearranged binary data, and
        corresponding stateframedef table as text file.
           sflog = file name of stateframe log to read
           sfxml = file name of corresponding XML file
           outbinfile = file name of output binary data file
           outtabfile = file name of output table text file
    '''
    if sfxml:
        sf, version = rxml.xml_ptrs(sfxml)
    else:
        sf = None
        version = 0.0

    if sflog:
        try:
            f = open(sflog,'rb')
        except:
            print 'Could not open file',sflog,'-- Exiting.'
            return
    
        # Get binary size and check version number
        data = f.read(100)
        if stateframe.extract(data,['d',8]) != version:
            print 'Stateframe log file version does not match XML version. -- Exiting'
            return
        recsize = stateframe.extract(data,sf['Binsize'])
        f.close()
    else:
        # No log file specified, so we will try to read directly from ACC once per second
        # Read one as a test and get its version number
        # Read from ACC
        accini = stateframe.rd_ACCfile()
        data, msg = stateframe.get_stateframe(accini)
        version = stateframe.extract(data,['d',8])
        

    # Parse the stateframe dictionary and generate the brange and outlist dicts
    brange, outlist = sfdef(sf)
    # Calculate the startbytes in the list -- modifies outlist in place
    startbyte(outlist)

    stdout = sys.stdout  # Save current stdout
    if outtabfile:
        # Write the table info to the given file name -- just sets stdout to the file,
        # writes it, and resets stdout
        try:
            sys.stdout = open(outtabfile,'w')
        except:
            print 'Could not redirect STDOUT to',outtabfile,' -- Will print to screen'
            sys.stdout = stdout

    outlist2table(outlist,version)
    sys.stdout = stdout   # Reset to standard stdout

    if outbinfile:
        try:
            g = open(outbinfile,'wb')
        except:
            print 'Could not open file',outbinfile,'for writing. -- Exiting.'
            return
        if sflog:
            # Read from log file
            f = open(sflog,'rb')
            while 1:
                # Read and rearrange 1000 records
                try:
                    indata = f.read(recsize)
                    outdata = transmogrify(indata,brange)
                    g.write(outdata)
                except:
                    f.close()
                    g.close()
                    return
        else:
            # Read from ACC
            accini = stateframe.rd_ACCfile()
            for i in range(60):
                # Read from ACC and rearrange 60 records -- takes 1 min
                indata, msg = stateframe.get_stateframe(accini)
                outdata = transmogrify(indata,brange)
                g.write(outdata)
                time.sleep(1)
            g.close()

    return            

#=============== load_deftable ===============
def load_deftable(xml_file=None,sdict=None,version=None):
    ''' Loads the named xml_file (or equivalent descriptive dictionary) 
        into the stateframe or scan header definition table, depending on 
        the type of xml file or dictionary provided.  Returns True if success, 
        or False if an error. It reads the current definition table to see 
        if this definition is already in there, and bails with a warning if so
    '''
    if xml_file and (not os.path.isfile(xml_file)):
        print 'Error: Named xml file',xml_file,'not found.'
        return False
    elif xml_file:
        try:
            sfname = os.path.basename(xml_file).split('.')[0][:10] == 'stateframe'
        except:
            sfname = False
        try:
            shname = os.path.basename(xml_file).split('.')[0][:11] == 'scan_header'
        except:
            shname = False
        if sfname or shname:
            sdict, version = rxml.xml_ptrs(xml_file)
    if sdict and version:
        brange, outlist = sfdef(sdict)
        startbyte(outlist)
        tbl = outlist2table(outlist,version)
        # Everything worked so far, so connect to database
        cnxn = pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;")
        cursor = cnxn.cursor()
        if 'Schedule' in sdict:
            # Case of a stateframe
            tblname = 'StateFrameDef'
        else:
            # Case of a scanheader
            tblname = 'ScanHeaderDef'
        # Check if this version is already entered
        print "select * from " + tblname + " where Version=" + str(int(version))
        cursor.execute("select * from " + tblname + " where Version=" + str(int(version)))
        rows = cursor.fetchall()
        if len(rows) == 0:
            for params in tbl:
                cursor.execute("insert into " + tblname + " (Status, Version, Dimension, "
                                  + "DataType, FieldBytes, DimOffset, StartByte, FieldNum, " 
                                  + "FieldName) values ( 0, ?, ?, ?, ?, ?, ?, ?, ?)", params)
            cursor.execute("update "+tblname+" set status=1 where Version='" + str(int(version)) + "'")
            cnxn.commit()
        else:
            print 'Warning: The definition for version',version,'already exists in',tblname
        cursor.close()
        del cursor
        cnxn.close()
        return True
    else:
        print 'Error: Bad sdict=',sdict,'or version=',version
        return False
    
def drop_deftable(version):
    ''' Drops ALL traces of a given version of a stateframe definition.
        Use with CAUTION!!!
        
        Requests confirmation from the keyboard.
    '''
    import dbutil as db
    cursor = db.get_cursor()
    # Get all table and view names from this version
    query = 'select * from information_schema.tables where table_name like "fV'+str(int(version))+'%"'
    result, msg = db.do_query(cursor, query)
    if msg == 'Success':
        names = result['TABLE_NAME']
        print 'You are about to permanently delete the following tables:'
        for name in names:
            print '    ',name
        ans = raw_input('Are you sure? [y/n]: ').upper()
        if ans != 'Y':
            print 'Action canceled by user'
            return False
        # Delete the version from the stateframedef table
        query = 'delete from StateFrameDef where Version='+str(int(version))
        r, msg = db.do_query(cursor, query)
        if msg != 'Success':
            print 'Error, could not delete from stateframedef for version:',version
            print msg
            return False
        # Loop over table and view names, dropping each in turn
        for name in names:
            query = 'drop table '+name
            r, msg = db.do_query(cursor, query)
            if msg != 'Success':
                print 'Error dropping table',name
                print msg
        # Drop Bin Reader procedure
        query = 'drop proc ov_fBinReader_V'+str(int(version))
        r, msg = db.do_query(cursor, query)
        if msg != 'Success':
            print 'Error, could not delete Bin Reader procedure for version:',version
            print msg
            return False
    else:
        return False
    print 'Successfully dropped all existing tables and table definition for version', version
    return True
        
#=============== reload_deftables ===============
def reload_deftables(tbldir='/data/eovsa/stateframe_logs'):
    ''' Clears the stateframe definition tables and reloads them from the defining xml tables
        NB: This only works on helios!
    '''
    if tbldir is None:
        print 'Error: no directory for xml tables was provided.'
#    files = glob.glob(os.path.join(tbldir,'*00.xml')
    logfiles = glob.glob(os.path.join(tbldir,'*0.log'))
    # Update tables only for log file versions that exist
    loglist = []
    for file in logfiles:
        loglist.append(file.split('_v')[1].split('.')[0])
    uniq = list(set(loglist))
    uniq.sort()
    for ver in uniq:
        xml_file = glob.glob(os.path.join(tbldir,'stateframe_v'+ver+'.00.xml'))
        if xml_file != []:
            load_deftable(xml_file[0])
    
#=============== log2sql ===============
def log2sql(log_file=None):
    ''' Transfers the named stateframe log file to the SQL database.  This transfer can
        take a long time, so this should allow interruption of the transfer, and then
        a subsequent call on the same log file should find the place where it left off to
        resume the transfer.
    '''
    from util import Time
    import traceback
    
    if log_file is None:
        print 'Error: a stateframe log filename must be provided.'
        return False
    if not os.path.isfile(log_file):
        print 'Error: Named stateframe log file',log_file,'not found.'
        return False
    # Log file basename is expected to be in format 'sf_yyyymmdd_vxx.0.log', 
    # where xx is the version number
    basename = os.path.basename(log_file) 
    logname = basename.split('_')
    if logname[0] == 'sf':
        sfdate = logname[1]
        try:
            sftime = datetime.datetime.strptime(logname[1],'%Y%m%d')
            t = Time(str(sftime))
            sftimestamp = int(t.lv + 0.5)  # Start timestamp, as nearest integer
            sfver = int(logname[2].split('.')[0][1:])
        except:
            print 'Error: File ',basename,'does not have expected basename format sf_yyyymmdd_vxx.0.log'
            return False
    else:
        return False
    
    # At this point, the log file exists and the name is of the right format
    # Connect to the database and see if there are any data already for this date, and if so
    # determine the time range.
    with pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;") as cnxn:
        cursor = cnxn.cursor()
        tblname = 'fV'+str(sfver)+'_vD1'
        cursor.execute("select top 1 Timestamp from "+tblname+" where Timestamp between "+str(sftimestamp)+" and "+str(sftimestamp+86400-2)+" order by Timestamp desc")
        rows = cursor.fetchall()
        if len(rows) == 1:
            # There are data for this date, and rows[1].Timestamp should be the last time entry,
            # so start at last time entry + 1 s
            try:
                sftimestamp2 = int(rows[0].Timestamp + 1)
            except:
                print 'Error: Unexpected data returned from database.  Returned value:',rows[0]
                return False
        elif len(rows) > 1:
            print 'Error: Unexpected data returned from database.'
            return False
        else:
            # No data returned from call, which means we should start with current value of sftimestamp
            pass
    
        # We now know where to start, so open log file and read to start of data
        # First need to find out record length
        f = open(log_file,'rb')
        buf = f.read(32)
        recsize = struct.unpack_from('i', buf, 16)[0]
        version = struct.unpack_from('d', buf, 8)[0]
        f.close()
        if int(version) != sfver:
            print 'Error: Version in file name is',sfver,'but version in file itself is',int(version)
            return False
            
        # We need the "brange" variable, which is used by transmogrify() to reformat the binary data.
        # Therefore, the defining stateframe XML file is needed.        
        # The correct XML file for this version must exist in the same directory as the log file
        xml_file = os.path.dirname(log_file)+'/'+'stateframe_v'+str(sfver)+'.00.xml'
        if not os.path.isfile(xml_file):
            print 'Error: Stateframe xml file',xml_file,'not found.'
            return False        
        sf, version = rxml.xml_ptrs(xml_file)
        brange, outlist = sfdef(sf)
        lineno = 0
        with open(log_file,'rb') as f:
            bufin = f.read(recsize)
            while len(bufin) == recsize:
                lineno += 1
                if struct.unpack_from('d', bufin, 0)[0] >= sftimestamp:
                    # This is new data, so write to database
                    bufout = transmogrify(bufin, brange)
                    try:
                        cursor.execute('insert into fBin (Bin) values (?)',pyodbc.Binary(bufout))
                        print 'Record '+str(lineno)+' successfully written\r',
                        cnxn.commit()
                    except Exception:
                        # An exception could be an error, or just that the entry was already inserted
                        traceback.print_exc()
                        pass
                bufin = f.read(recsize)
    print '\n'
    return True

#=============== log2sql ===============
def badlog2sql(log_file=None):
    ''' This version checks for bad (short) records in the log file (caused by a bug
        in the schedule in March/April 2022) and skips any short records.
        
        Transfers the named stateframe log file to the SQL database.  This transfer can
        take a long time, so this should allow interruption of the transfer, and then
        a subsequent call on the same log file should find the place where it left off to
        resume the transfer.
    '''
    from util import Time
    import traceback
    
    if log_file is None:
        print 'Error: a stateframe log filename must be provided.'
        return False
    if not os.path.isfile(log_file):
        print 'Error: Named stateframe log file',log_file,'not found.'
        return False
    # Log file basename is expected to be in format 'sf_yyyymmdd_vxx.0.log', 
    # where xx is the version number
    basename = os.path.basename(log_file) 
    logname = basename.split('_')
    if logname[0] == 'sf':
        sfdate = logname[1]
        try:
            sftime = datetime.datetime.strptime(logname[1],'%Y%m%d')
            t = Time(str(sftime))
            sftimestamp = int(t.lv + 0.5)  # Start timestamp, as nearest integer
            sfver = int(logname[2].split('.')[0][1:])
        except:
            print 'Error: File ',basename,'does not have expected basename format sf_yyyymmdd_vxx.0.log'
            return False
    else:
        return False
    
    # At this point, the log file exists and the name is of the right format
    # Connect to the database and see if there are any data already for this date, and if so
    # determine the time range.
    with pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;") as cnxn:
        cursor = cnxn.cursor()
        tblname = 'fV'+str(sfver)+'_vD1'
        cursor.execute("select top 1 Timestamp from "+tblname+" where Timestamp between "+str(sftimestamp)+" and "+str(sftimestamp+86400-2)+" order by Timestamp desc")
        rows = cursor.fetchall()
        if len(rows) == 1:
            # There are data for this date, and rows[1].Timestamp should be the last time entry,
            # so start at last time entry + 1 s
            try:
                sftimestamp2 = int(rows[0].Timestamp + 1)
            except:
                print 'Error: Unexpected data returned from database.  Returned value:',rows[0]
                return False
        elif len(rows) > 1:
            print 'Error: Unexpected data returned from database.'
            return False
        else:
            # No data returned from call, which means we should start with current value of sftimestamp
            pass
    
        # We now know where to start, so open log file and read to start of data
        # First need to find out record length
        f = open(log_file,'rb')
        buf = f.read(32)
        recsize = struct.unpack_from('i', buf, 16)[0]
        version = struct.unpack_from('d', buf, 8)[0]
        srchstr = buf[8:24]   # Binary-coded version string + recsize
        f.close()
        if int(version) != sfver:
            print 'Error: Version in file name is',sfver,'but version in file itself is',int(version)
            return False
            
        # We need the "brange" variable, which is used by transmogrify() to reformat the binary data.
        # Therefore, the defining stateframe XML file is needed.        
        # The correct XML file for this version must exist in the same directory as the log file
        xml_file = os.path.dirname(log_file)+'/'+'stateframe_v'+str(sfver)+'.00.xml'
        if not os.path.isfile(xml_file):
            print 'Error: Stateframe xml file',xml_file,'not found.'
            return False        
        sf, version = rxml.xml_ptrs(xml_file)
        brange, outlist = sfdef(sf)
        lineno = 0
        fptr = 0
        with open(log_file,'rb') as f:
            bufin = f.read(recsize)
            while len(bufin) == recsize:
                lineno += 1
                fptr += recsize   # Current byte in file
                pos = bufin[16:].find(srchstr) + 16    # See if version+recsize occurs elsewhere in buffer => short buffer
                if pos > 15:
                    # Bad record, so seek to that position-8 (back up)
                    print 'Bad rec at', lineno, 'skipped'
                    fptr -= recsize - pos + 8
                    f.seek(fptr)
                else:
                    if struct.unpack_from('d', bufin, 0)[0] >= sftimestamp:
                        # This is new data, so write to database
                        bufout = transmogrify(bufin, brange)
                        try:
                            cursor.execute('insert into fBin (Bin) values (?)',pyodbc.Binary(bufout))
                            print 'Record '+str(lineno)+' successfully written\r',
                            cnxn.commit()
                        except Exception:
                            # An exception could be an error, or just that the entry was already inserted
                            traceback.print_exc()
                            pass
                bufin = f.read(recsize)
    print '\n'
    return True

def scanheader_log2sql(log_file):
    from util import Time
    import traceback
    
    if log_file is None:
        print 'Error: a scanheader log filename must be provided.'
        return False
    if not os.path.isfile(log_file):
        print 'Error: Named scanheader log file',log_file,'not found.'
        return False
    # Log file basename is expected to be in format 'sh_yyyymmdd_vxx.0.log', 
    # where xx is the version number
    basename = os.path.basename(log_file) 
    logname = basename.split('_')
    if logname[0] == 'sh':
        sfdate = logname[1]
        try:
            shtime = datetime.datetime.strptime(logname[1],'%Y%m%d')
            t = Time(str(shtime))
            shtimestamp = int(t.lv + 0.5)  # Start timestamp, as nearest integer
            shver = int(logname[2].split('.')[0][1:])
        except:
            print 'Error: File ',basename,'does not have expected basename format sh_yyyymmdd_vxx.0.log'
            return False
    else:
        return False
    
    # At this point, the log file exists and the name is of the right format
    # Connect to the database and see if there are any data already for this date, and if so
    # determine the time range.
    with pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;") as cnxn:
        cursor = cnxn.cursor()
        tblname = 'hV'+str(shver)+'_vD1'
        cursor.execute("select top 1 Timestamp from "+tblname+" where Timestamp between "+str(shtimestamp)+" and "+str(shtimestamp+86400-2)+" order by Timestamp desc")
        rows = cursor.fetchall()
        if len(rows) == 1:
            # There are data for this date, and rows[1].Timestamp should be the last time entry,
            # so start at last time entry + 1 s
            try:
                shtimestamp2 = int(rows[0].Timestamp + 1)
            except:
                print 'Error: Unexpected data returned from database.  Returned value:',rows[0]
                return False
        elif len(rows) > 1:
            print 'Error: Unexpected data returned from database.'
            return False
        else:
            # No data returned from call, which means we should start with current value of sftimestamp
            pass
    
        # We now know where to start, so open sample data file and get its length
        template_file = os.path.dirname(log_file)+'/'+'scan_header.dat'
        # First need to find out record length
        f = open(template_file,'rb')
        buf = f.read()
        recsize = len(buf)
        version = struct.unpack_from('d', buf, 8)[0]
        srchstr = buf[8:16]   # Binary-coded version string
        f.close()
        if int(version) != shver:
            print 'Error: Version in file name is',shver,'but version in file itself is',int(version)
            return False
            
        # We need the "brange" variable, which is used by transmogrify() to reformat the binary data.
        # Therefore, the defining stateframe XML file is needed.        
        # The correct XML file for this version must exist in the same directory as the log file
        xml_file = os.path.dirname(log_file)+'/'+'scanheader_v'+str(shver)+'.00.xml'
        if not os.path.isfile(xml_file):
            print 'Error: Scanheader xml file',xml_file,'not found.'
            return False
        sh, version = rxml.xml_ptrs(xml_file)
        brange, outlist = sfdef(sh)
        lineno = 0
        fptr = 0
#        import pdb; pdb.set_trace()
        with open(log_file,'rb') as f:
            bufin = f.read(recsize)
            while len(bufin) == recsize:
                lineno += 1
                fptr += recsize   # Current byte in file
                pos = bufin[16:].find(srchstr) + 16    # See if version+recsize occurs elsewhere in buffer => short buffer
                if pos > 15:
                    # Bad record, so seek to that position-8 (back up)
                    print 'Bad rec at', lineno, 'truncated'
                    fptr -= recsize - pos + 8
                    f.seek(fptr)
                # else:
                if struct.unpack_from('d', bufin, 8)[0] >= version:
                    # This is new data, so write to database
                    bufout = transmogrify(bufin, brange)
                    try:
                        cursor.execute('insert into hBin (Bin) values (?)',pyodbc.Binary(bufout))
                        print 'Record '+str(lineno)+' successfully written\r',
                        cnxn.commit()
                    except Exception:
                        # An exception could be an error, or just that the entry was already inserted
                        traceback.print_exc()
                        pass
                else:
                    print 'Error in record',lineno,'so not sent.'
                bufin = f.read(recsize)
    print '\n'
    return True

#=============== acc2sql ===============
def acc2sql():
    ''' This is just a test version to read the stateframe once a second from
        the ACC and send it to the SQL server.  A more complete version of this
        is implemented in schedule.py for "production" work.
    '''
    # Get stateframe structure and version
    accini = stateframe.rd_ACCfile()
    sf_version = accini['version']
    brange, outlist = sfdef(accini['sf'])
    with pyodbc.connect("DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                             DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;") as cnxn:
        cursor = cnxn.cursor()
        lineno = 0
        while 1:
            lineno += 1
            data, msg = stateframe.get_stateframe(accini)
            version = stateframe.extract(data,['d',8])
            if version == sf_version:
                bufout = transmogrify(data, brange)
                try:
                    cursor.execute('insert into fBin (Bin) values (?)',pyodbc.Binary(bufout))
                    print 'Record '+str(lineno)+' successfully written\r',
                    cnxn.commit()
                except:
                    # An exception could be an error, or just that the entry was already inserted
                    pass
            else:
                print 'Error: Incompatible version in stateframe.'
                break
            time.sleep(1)

