#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#

#from lxml import etree
import xml.etree.ElementTree as etree
import numpy as np
import copy
import struct
import urllib2
import socket

def handle_cluster(child):
    '''This element of the XML tree is the head of a Cluster.  Step through
       each element of the branch and return the keys, the empty dictionary, 
       and the fmt string.  This routine is reentrant'''
    # Clusters have a name, a NumElts, and one or more objects
    c = list(child)
    fmt = ''
    if c[0].tag == "Name":
        keys = [c[0].text]
        mydict = {}    # Start an empty dictionary
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None
    if c[0].tag == "NumElts":
        n = int(c[0].text)
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None
    # Loop through all items in the cluster
    for i in range(n):
        datatype = c[0].tag
        if datatype == "Array":
            ch = c[0]
            key, arr, dims, fmt1 = handle_array(ch)
            keys += key
            fmt += fmt1
            mydict.update({key[0]:arr})  # Add a list to the dictionary
            # Have to figure out what to do with dims -- maybe ignore and use fmt code?
            c.pop(0)
        elif datatype == "Cluster":
            ch = c[0]
            key, newdict, fmt1 = handle_cluster(ch)
            if not key is None:
                keys += key
            fmt += fmt1
            mydict.update({key[0]:newdict})  # Add a dictionary to the dictionary
            c.pop(0)
        else:
            ch = c[0]
            key, fmt1 = handle_item(ch)
            keys += key
            fmt += fmt1
            mydict.update({key[0]:0})   # Add an item to the dictionary
            c.pop(0)
    return keys, mydict, fmt

def handle_array(child):
    '''This element of the XML tree is an Array.  Step through the items of the
       array (which may contain clusters, other arrays, etc.) and return the keys,
       array, dimensions of the array, and fmt string.  This routine is reentrant.
    '''
    # Arrays have a name and one or more dimension statements, then one or more objects
    c = list(child)
    if c[0].tag == "Name":
        keys = [c[0].text]
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None, None, None
    # Handle up to four levels of dimension
    d1, d2, d3, d4 = 1, 1, 1, 1
    if c[0].tag == "Dimsize":
        d1 = int(c[0].text)
        fmt = 'I'
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None
    if c[0].tag == "Dimsize":
        d2 = int(c[0].text)
        fmt += 'I'
        c.pop(0)
    if d2 != 1:
        if c[0].tag == "Dimsize":
            d3 = int(c[0].text)
            fmt += 'I'
            c.pop(0)
        if d3 != 1:
            if c[0].tag == "Dimsize":
                d4 = int(c[0].text)
                fmt += 'I'
                c.pop(0)
    dims = [d1, d2, d3, d4]
    datatype = c[0].tag
    dtype_dict = {'U8':'s','B8':'B','U16':'H','U32':'I','I16':'h','I32':'i','SGL':'f','DBL':'d'}
    fmt += str(d1*d2*d3*d4)+dtype_dict.get(datatype,'[')
    if datatype == "Cluster":
        ch = list(c[0])
        key, mydict, fmt1 = handle_cluster(ch)
        keys += key
        fmt += fmt1+']'
        arr = [mydict]  # Return cluster dictionary as 1-element list place holder
        c.pop(0)
    else:
        arr = dims   # Return list of dims as place holder
    return keys, arr, dims, fmt

def handle_item(c):
    '''This element of the XML tree is a simple, single item.  Simply return its key and fmt string.
    '''
    dtype_dict = {'U8':'s','B8':'B','U16':'H','U32':'I','I16':'h','I32':'i','SGL':'f','DBL':'d'}
    fmt = dtype_dict.get(c.tag,'*')
    if c[0].tag == "Name":
        key = [c[0].text]
    else:
        print 'Illegal format for item',c
        return None, None
    return key, fmt

def xml_read(filename=None):
    '''This is a pre-processing step for reading an XML file and creating a
       structure describing its contents.  It "walks" the XML tree and 
       recursively applies the routines
           handle_item()
           handle_array()
           handle_cluster()
       From an XML description of the stateframe it returns a dictionary with 
       key:value pairs or lists thereof, and a format string for decoding the 
       stateframe.  Items returned: 
           keys    list of strings that are variable names in the XML file,
                     in the order they appear in the XML file, to get around 
                     the fact that dictionary keys can appear in any order.
           mydict  a dictionary whose keys are now in random order, with values
                     that are place-holders for information that will be read
           fmt     a pseudo-Python struct string that contains non-Python 
                     '[' and ']' to indicate the extent of arrays of Clusters. 
    '''
    userpass = 'admin:observer@'

    if filename is None:
        f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/ni-rt/startup/stateframe.xml')
        if socket.gethostname() == 'helios':
            # If this is the OVSA machine, make a disk copy of stateframe.xml in the
            # current (dropbox) directory.  This will be used by other instances of
            # sf_display() on other machines that do not have access to acc.solar.pvt.
            lines = f.readlines()
            o = open('stateframe.xml','w')
            for line in lines:
                o.write(line+'\n')
            o.close()
            f.close()
            f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/ni-rt/startup/stateframe.xml')
    else: 
        f = open(filename)
    tree = etree.parse(f)
    f.close()

    root = tree.getroot()
    try:
        # Stateframe version number is supposed to be included in the second
        # element of the stateframe cluster (after the timestamp)
        version = float(root[3][1].text)
    except:
        # It seems not to be there, so set version to 3.0
        version = 3.0

    keys, mydict, fmt = handle_cluster(root)
    return keys, mydict, fmt, version

def handle_key(keys, dictlist, fmt, off):
    '''This routine steps through the keys list from the preprocessing step
       and creates dictionary entries for each data element (be it an array,
       a cluster, or a single value).  Each cluster is represented by a
       separate dictionary, and nested dictionaries are listed in dictlist,
       with the currently active, inner-most dictionary being the last in 
       the list.
       For example, the weather value Temperature is in the 3rd-level 
       dictionary 
          sf['Schedule']['Data']['Weather']['Temperature']. 
       For this case, dictlist would be [{Schedule},{Data},{Weather}], and 
       dictlist.pop() would present the innermost {Weather} dictionary for 
       processing.  
       Each pass through this recursive routine returns its input arguments,
       with keys and fmt pared down, and off increased, as each key is 
       handled, until fmt is returned empty.  The final dictlist is a single 
       dictionary with its values being lists of individual fmt strings and 
       offsets into the data buffer.
    '''
    key = keys.pop(0)
    if key is None:
        #Skip any "None" key
        return keys, dictlist, fmt, off
    mydict = dictlist.pop()      # Get the inner-most dict for processing
    try:
        val = mydict[key]
    except:
        # We must be done with this dictionary, so go back and try again
        # Note that dictlist is one item shorter than before
        keys = [key] + keys  # Put key back in list
        return keys, dictlist, fmt, off
    valtype = type(val)
    if valtype == int or valtype == float:
        # This is just a single value
        f = fmt[0]
        fmt = fmt[len(f):]
        try:
            mydict[key] = [f, off]   # Assign fmt, offset pair as value to key
            # Increment off by number of bytes taken by f            
            off += struct.calcsize(f)
        except:
            print key, fmt
        dictlist.append(mydict)   # Put original dictionary back
    elif valtype == list:
        # This is an array of values
        if type(val[0]) == int:
            # This is a simple array of numbers
            dims = val  # List of dimensions from XML file
            arrsiz = 1  # Total array size (product of dimensions)
            ndim = 0    # Number of dimensions 
            for dim in dims:
               arrsiz *= dim
               if dim != 1:
                   # Count only non-unity dimensions
                   ndim += 1
            dims = dims[:ndim]
            # Read dimensions of array
            while fmt[0] == 'I':
                f = fmt[0]
                fmt = fmt[len(f):]
                # Skip dimension variables
                # Increment off by number of bytes taken by f
                off += struct.calcsize(f)
            f = fmt[:len(str(arrsiz))+1]
            fmt = fmt[len(f):]
            if f[-1] == 's':
                mydict[key] = [f, off]  # If f is 'string', do not save dims
            else:
                mydict[key] = [f, off, dims]   # Assign fmt, offset and dims as value to key
            # Increment off by number of bytes taken by f (to prepare for next iteration)
            off += struct.calcsize(f)
            dictlist.append(mydict)    # Put original dictionary back
        else:
            # This is something more complicated (e.g. an array of dicts)
            if type(val[0]) == dict:
#                dims = []   # List of dimensions
                arrsiz = 1  # Total array size (product of dimensions)
                dictarr = []
                # Read dimensions of array
                while fmt[0] == 'I':
                    f = fmt[0]
                    fmt = fmt[len(f):]
                    # Skip dimension variables
                    # Increment off by number of bytes taken by f
                    off += struct.calcsize(f)
#                    dims.append(vals[0])
#                    arrsiz *= vals[0]
                # Extract array size (number just before '[')
                arrsiz = int(fmt[:fmt.index('[')])   # This may break for some "valid" XML files
                # Extract format string from between [] brackets, which will be applied repeatedly
                newfmt = fmt[fmt.index('[')+1:fmt.index(']')]
                fmt = fmt[fmt.index(']')+1:]  # Remaining fmt after closing ] bracket
                # Loop over total number of dicts in array
                for j in range(arrsiz):
                    newdictarr = [copy.deepcopy(val[0])]
                    tmpfmt = newfmt
                    tmpkeys = copy.deepcopy(keys)
                    while tmpfmt != '':
                        tmpkeys, newdictarr, tmpfmt, off = handle_key(tmpkeys, newdictarr, tmpfmt, off)
                    while len(newdictarr) > 1:
                        newdictarr.pop()  # Remove all but the original dictionary
                    # Dictionary is all filled in, and only one copy remains.  Copy it to dictarray being assembled
                    dictarr.append(copy.deepcopy(newdictarr.pop()))
                keys = tmpkeys
            mydict[key] = dictarr    # Assign array of dicts to mydict key
            dictlist.append(mydict)  # Put original dictionary back
    elif valtype == dict:
        # This is a dictionary.
        dictlist.append(mydict)  # Put original dictionary back
        dictlist.append(val)     # Add new dictionary
    else:
        print 'Unknown value type',valtype
    return keys, dictlist, fmt, off

def xml_ptrs(filename=None):
    '''Reads XML file filename, and creates a nested dictionary structure with
       keys from the XML file and values that are lists of the form 
       [fmt, offset], where fmt is the struct format string for a data element
       and offset is the byte offset into the stateframe data buffer of the
       start of the data element.  An example of use is the AzimuthError1
       element from antenna 1, which would be read from a stateframe log file
       using the following code:
           f = open('stateframe.log')
           sf, version = xml_ptrs('stateframe.xml')
           data = f.read(reclen)   # Read the data buffer from file-handle f
           fmt, off = sf['Antenna'][0]['Controller']['AzimuthError1']
           azerr = struct.unpack_from(fmt,data,off)
       Also returned is the version variable, which is the currently read XML file
       version for comparison with the version number in the binary data.
    '''
    inkeys, indict, infmt, version = xml_read(filename)   # Pre-processing step
    keys = copy.deepcopy(inkeys)
    mydict = copy.deepcopy(indict)
    fmt = infmt
    keys.pop(0)
    dictlist = [mydict]  # A list of dictionaries.  The one at the end is the one currently being manipulated
    off = 0
    while fmt != '':
        keys, dictlist, fmt, off = handle_key(keys, dictlist, fmt, off)
    mydict = dictlist.pop()   # Should be the original, but now updated dictionary
    if dictlist != []:
        mydict = dictlist.pop()   # Should be the original, but now updated dictionary
    return mydict, version

