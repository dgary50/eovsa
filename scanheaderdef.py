import read_xml2, stateframe, struct, numpy

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
                    # that there are two special-case exceptions to
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
    
def rd_sfdef(sh=None):
    ''' Create initial list of dictionaries representing each variable
        in the given scan header dictionary sh.  If not given, sh is
        obtained from the current version.
    '''
    if sh is None:
        sh, version = read_xml2.xml_ptrs('/tmp/scan_header.xml')

    mylist = walk_keys(sh,'Stateframe')
    return mylist

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
        into lists of a given dimension
    '''
    dim = mylist[0]['dimension']
    lists = []
    sublist = []
    for d in mylist:
        if d['dimension'] == dim:
            sublist += [d]
        else:
            lists += [sublist]
            sublist = [d]
            dim = d['dimension']
    lists += [sublist]            
    return lists

def sfdef(sf=None):
    ''' Main routine to generate a dictionary list representing the
        stateframedef table (outlist), and a dictionary list of
        rearrangements to the stateframe binary data, brange.

        Variable sf is an optional stateframe dictionary.  If omitted, the
        current stateframe.xml file is used generate it.
    '''
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
            if len(item['pytype']) > 1:
                # Case of "in-place" array
                ebyte = item['startbyte']+item['fieldbytes']*item['dimension']
            else:
                ebyte = item['startbyte']+item['fieldbytes']
            if ebyte != slist[i+1]['startbyte']:
                # There is a gap between this and next item, but this might
                # be due to dimensionality of this section.  Add bytes for
                # dimensionality and check again.
                if len(item['pytype']) > 1:
                    # Case of "in-place" array -- do nothing
                    pass
                else:
                    ebyte = (ebyte - sbyte)*item['dimension'] + sbyte
                if ebyte != slist[i+1]['startbyte']:
                    # This is indeed a gap, so save byte range and start a
                    # new byte range
                    brange.append({'sbyte':sbyte,'ebyte':ebyte})
                    sbyte = slist[i+1]['startbyte']
        ebyte = slist[-1]['startbyte']+slist[-1]['fieldbytes']
        # Take care of dimensionality of last item
        ebyte = (ebyte - sbyte)*slist[-1]['dimension'] + sbyte
        brange.append({'sbyte':sbyte,'ebyte':ebyte})
    # At this point, we should have a complete map of the ordered byte
    # ranges into original binary data, in variable brange.
    return brange, outlist

def startbyte(outlist):
    # The input startbyte has to be transferred to origbyte, and
    # we have to do another pass through to calculate output startbyte
    curpos = 0
    dim = outlist[0]['dimension']
    offset = outlist[0]['dimoffset']
    for i in range(len(outlist)):
        outlist[i].update({'origbyte':outlist[i]['startbyte']})
        outlist[i].update({'startbyte':curpos})
        if outlist[i]['dimension'] == 1:
            # For dimension 1 items, dimoffset is irrelevant, so set to 0
            outlist[i]['dimoffset'] = 0
        # If the current line has a dimension change or a dimoffset change...
        if (outlist[i]['dimension'] != dim) or (outlist[i]['dimoffset'] != offset) or (len(outlist[i]['pytype']) > 1):
            # Update the current position by n-1 x offset of previous line
            curpos += outlist[i-1]['dimoffset']*(dim-1)
            outlist[i].update({'startbyte':curpos})
            dim = outlist[i]['dimension']
            offset = outlist[i]['dimoffset']
        curpos += outlist[i]['fieldbytes']

def outlist2table(outlist,version):
    ''' Convert the completed outlist to an SQL table, which is the actual
        stateframedef table.  For now, this just prints a formatted table,
        but it will eventually send SQL commands to create the table.
    '''
    idx = 1
    dim = 1
    # If version is an integer (i.e. no float remainder) then write it out as
    # an integer (this converts, e.g., 25.0 to 25, but leaves 25.1 as is)
    ver = version
    if (version % 1) == 0:
        ver = int(version)
    for i,item in enumerate(outlist):
        if item['dimension'] != dim:
            # This generates the "fieldnum," which is just the ordinal position in a table
            # of a given dimension.  It resets to 1 whenever a dimension change occurs in
            # the ordered table.
            idx = 1
            dim = item['dimension']
        print '{:8} {:8} {:8} {:8} {:8} {:3} {:8} {:8} {:8}    {:}'.format(i+1,0,ver,item['dimension'],item['datatype'],item['fieldbytes'],item['dimoffset'],item['startbyte']+1,idx,item['fieldname'])
        idx += 1
        
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

import read_xml2 as rxml
import sys

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
            for i in range(1000):
                # Read and rearrange 1000 records
                try:
                    indata = f.read(recsize)
                    outdata = transmogrify(indata,brange)
                    g.write(outdata)
                except:
                    break
            f.close()
            g.close()
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
