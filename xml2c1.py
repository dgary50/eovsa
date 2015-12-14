#
# History:
#   2015-Oct-09  DG
#     In c_handle_array(), fixed case where arrays of clusters were repeated, but 
#     not given unique names
#
from lxml import etree
import numpy as np
import copy
import struct
import urllib2
import time

global indent
indent = 0

# slist is for structures that have been defined
global slist
slist = ['dude']


def c_handle_cluster(child):
    '''This element of the XML tree is the head of a Cluster.  Step through
       each element of the branch and return the keys, the empty dictionary, 
       and the fmt string.  This routine is reentrant'''
    global indent
    global slist
    # Clusters have a name, a NumElts, and one or more objects
    c = list(child)
    if c[0].tag == "Name":
        sname = c[0].text  # Name of structure
        lines = []
        if sname:
            #here find if sname has been used before, if so then add a '1' to the name to make it unique
            testname = 'x'+sname.lower()
            while slist.count(testname) > 0:
                testname = testname+'1'
            lines = ['struct '+testname+' {'] # Start struct line list
            slist = slist+[testname]
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
        objecttype = c[0].tag
        if objecttype == "Array":
            ch = c[0]
            name, newlines, dims, datatype = c_handle_array(ch)
            # Create line declaring variable 'size_'+name for holding dimensions
            line = ' '*3+'int size_'+name
            if len(dims) > 1:
                line += '['+str(len(dims))+']'
            line += ';'
            lines.append(line)
            # Create line declaring array variable itself
            # If datatype is a structure, we will later insert the structure
            # definition between datatype and name
            line = ' '*3+datatype+' '
            if newlines:
                # If this is a structure, we put the declaration line
                # and then append the newlines, before inserting name
                lines.append(line)
                line = ' '*6+'} '
                for i in range(len(newlines)):
                    newlines[i] = ' '*3+newlines[i]
                lines += newlines
            # Insert the name and dimensions
            line += name+'['
            # Replace any commas in dim string with ']['
            for dim in dims:
                line += str(dim)+']['
            line = line[:len(line)-2]+'];'
            lines.append(line)
            c.pop(0)
        elif objecttype == "Cluster":
            ch = c[0]
            name, newlines = c_handle_cluster(ch)
            if name is None:
                name = arrayname
            # Add indentation spaces for nesting level
            for i in range(len(newlines)):
                newlines[i] = ' '*3+newlines[i]
            lines += newlines
            #Be sure that name has no spaces
            name1=str.split(name)
            if len(name1) > 1:
                name=''.join(name1)
            #Be sure that name has no dots
            name1=str.split(name, '.')
            if len(name1) > 1:
                name=''.join(name1)

            lines += [' '*3+'   } '+name+';']
            c.pop(0)
        else:
            ch = c[0]
            name, datatype = c_handle_item(ch)
            # Add a simple line entry with datatype followed by name of variable
            #Be sure that name has no spaces
            name1=str.split(name)
            if len(name1) > 1:
                name=''.join(name1)
            #Be sure that name has no dots
            name1=str.split(name, '.')
            if len(name1) > 1:
                name=''.join(name1)
            lines += [' '*3+datatype+' '+name+';']
            c.pop(0)
    return sname, lines

def c_handle_array(child):
    '''This element of the XML tree is an Array.  Step through the items of the
       array (which may contain clusters, other arrays, etc.) and return the keys,
       array, dimensions of the array, and fmt string.  This routine is reentrant.
    '''
    global slist
    # Arrays have a name and one or more dimension statements, then one or more objects
    c = list(child)
    if c[0].tag == "Name":
        name = c[0].text  # Name of Array
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None, None, None
    # Handle up to four levels of dimension
    dims = []
    if c[0].tag == "Dimsize":
        dims.append(int(c[0].text))
        c.pop(0)
    else:
        print 'Illegal format for item',child
        return None, None, None, None
    if c[0].tag == "Dimsize":
        dims.append(int(c[0].text))
        c.pop(0)
    if c[0].tag == "Dimsize":
        dims.append(int(c[0].text))
        c.pop(0)
    if c[0].tag == "Dimsize":
        dims.append(int(c[0].text))
        c.pop(0)

    # It turns out thet the dims need reversing, do that here, jmm, 14-sep-2013
    dims.reverse()

    objecttype = c[0].tag
    dtype_dict = {'U8':'char','B8':'unsigned char','U16':'unsigned short int','U32':'unsigned int',
                  'I16':'short int','I32':'int','SGL':'float','DBL':'double'}
    datatype = dtype_dict.get(objecttype,'*')
    if objecttype == "Cluster":
        ch = list(c[0])
        name2, lines = c_handle_cluster(ch)
        if name2:
            name = name2
        #here find if name has been used before, if so then add a '1' to the name to make it unique
        testname = 'x'+name.lower()
        while slist.count(testname) > 0:
            testname = testname+'1'
        datatype = 'struct '+testname+' {' # Start struct line list
        slist = slist+[testname]
        c.pop(0)
    else:
        lines = []
    return name, lines, dims, datatype

def c_handle_item(c):
    '''This element of the XML tree is a simple, single item.  Simply return its key and fmt string.
    '''
    dtype_dict = {'U8':'char','B8':'char','U16':'unsigned short int','U32':'unsigned int',
                  'I16':'short int','I32':'int','SGL':'float','DBL':'double'}
    datatype = dtype_dict.get(c.tag,'*')
    if c[0].tag == "Name":
        name = c[0].text
    else:
        print 'Illegal format for item',c
        return None, None
    return name, datatype

def c_xml_read_stateframe(filename=None):
    '''Read an XML description of the stateframe and return lines for output
    to a C-language header file.  This header file can be included in a
    compilation to allow reading the binary data corresponding to the XML file
    and memory mapping it to the structure.
    '''
    if filename is None:
        f = urllib2.urlopen('ftp://admin:observer@acc.solar.pvt/ni-rt/startup/stateframe.xml')
    else: 
        f = open(filename)
    tree = etree.parse(f)
    f.close()

    root = tree.getroot()
    name, lines = c_handle_cluster(root)
    lines.append('   } '+name+';')

    #grab version here too
    if filename is None:
        f = urllib2.urlopen('ftp://admin:observer@acc.solar.pvt/ni-rt/startup/stateframe.xml')
    else: 
        f = open(filename)
    
    p = f.readlines()
    f.close()

    #step through and find the version
    version = "NONE"
    j = -1
    while version == "NONE" and j < (len(p)-1):
        j = j+1
        tbp = str.find(p[j], 'Version')
        #If we have a version, then extract it
        if tbp != -1:
            j1 = j+1
            version_line = p[j1]
            #strip off <Val>, </Val>
            version = version_line[5:len(version_line)-8]
            
    
    return lines, version

def c_xml_read_scanheader(filename=None):
    '''Read an XML description of the scanheader and return lines for output
    to a C-language header file.  This header file can be included in a
    compilation to allow reading the binary data corresponding to the XML file
    and memory mapping it to the structure.
    '''
    if filename is None:
        f = urllib2.urlopen('ftp://admin:observer@acc.solar.pvt/parm/scan_header.xml')
    else: 
        f = open(filename)
    tree = etree.parse(f)
    f.close()

    root = tree.getroot()
    name, lines = c_handle_cluster(root)
    lines.append('   } '+name+';')

    #grab version here too
    if filename is None:
        f = urllib2.urlopen('ftp://admin:observer@acc.solar.pvt/parm/scan_header.xml')
    else: 
        f = open(filename)
    
    p = f.readlines()
    f.close()

    #step through and find the version
    version = "NONE"
    j = -1
    while version == "NONE" and j < (len(p)-1):
        j = j+1
        tbp = str.find(p[j], 'Version')
        #If we have a version, then extract it
        if tbp != -1:
            j1 = j+1
            version_line = p[j1]
            #strip off <Val>, </Val>
            version = version_line[5:len(version_line)-8]

    return lines, version


def c_struct_unnest(arrx1):
    '''Takes a C structure definintion for the stateframe, and unnests it into
    something that can be translated into a FORTRAN common block. Input is a 
    list of lines in the structure definition'''

    # create a 2-d list with each element
    nlines = len(arrx1)
    arrx0 = range(nlines)
    for j in range(nlines):
        arrx0[j] = str.split(arrx1[j], ' ')
    
    # str.split with a ' ' leaves null strings for spaces, we can
    # count to get indentation
    # for each line, find the first value and 
    # create a copy of arrx0 without the spaces in front
    level = range(nlines)
    first_x = range(nlines)
    for j in range(nlines):
        a = 0
        while arrx0[j][a] == '':
            a = a+1
        level[j] = a/3
        first_x[j] = arrx0[j][a]

    # find the start of the structures
    xxx = [-1]
    for j in range(nlines):
        if first_x[j] == 'struct':
            xxx = xxx+[j]

    xxx = xxx[1:]
        
    nstruct = len(xxx)
    xxx1 = range(nstruct)
    for j in range(nstruct):
        xxx1[j] = nlines-1

    nest_level = range(nstruct)
    for j in range(nstruct):
        j1 = xxx[j]
        nest_level[j] = level[j1]
        # here i need to find the next point where first_x is '}' and level is nest_level[j]+1
        while j1 < (nlines-1) and not(first_x[j1] == '}' and level[j1] == nest_level[j]+1):
            j1 = j1+1
        xxx1[j] = j1

    # Here create a list which for each line has either zero,
    # or for the positions xxx of the structure starts,
    # the value of the appropriate xxx1, the structure end
    stest = range(nlines)
    for j in range(nlines):
        stest[j] = 0

    for j in range(nstruct):
       for i in range(nlines):
           if i == xxx[j]:
               stest[i] = xxx1[j]

    # For each nest level, get the structure definitions first:
    nl = sorted(nest_level)
    max_nl = nl[len(nl)-1]+1 #Another random plus 1?
    # Need to reverse a list
    pppp = range(max_nl)
    qqqq = pppp[::-1]
    otp_string = [-1]
    for j in qqqq:
        n = [-1]
        for k in range(nstruct):
            if nest_level[k] == j:
                n = n+[k]
        
        n = n[1:]
        for k in range(len(n)):
            #I don't know why I need the plus 1 in here yet
            st_st = xxx[n[k]]
            en_st = xxx1[n[k]] 
            strdefk = arrx1[st_st:en_st+1]
            stestk = stest[st_st:en_st+1]
            #Unless nest_level is zero, strip the variable 
            #name from the last element
            if j > 0 :
                strdefk[len(strdefk)-1]='   };'

            # Build the output, strip out sub structure definitions
            nk = len(strdefk)
            in_substruct = 0
            end_substruct = nk-1
            for l in range(nk):
                #If i'm inside a substructure then do nothing
                if in_substruct == 0:
                    if l == 0:
                        temp_string = strdefk[l]
                        # remove leading spaces -- j is the nest level
                        otp_string = otp_string+[temp_string[3*j:]]        
                    elif l == nk-1:
                        temp_string = strdefk[l]
                        otp_string = otp_string+[temp_string]        
                    else:
                        if stestk[l] != 0:
                            # strip out the trailing bracket
                            temp_l = strdefk[l]
                            len_l = len(temp_l)
                            temp_l = temp_l[:len_l-2]
                            in_substruct = 1
                            end_substruct = stestk[l]-st_st
                        else:
                            temp_string = strdefk[l]
                            otp_string = otp_string+[temp_string[3*j:]]        
                else:
                    if l == end_substruct:
                        in_substruct = 0
                        temp_l1 = strdefk[l]
                        bp = str.find(temp_l1, '}')
                        temp_l1 = temp_l1[bp+1:]
                        temp_string = temp_l+temp_l1
                        otp_string = otp_string+[temp_string[3*j:]]        


# Done? Except that xattenuation is defined twice, 
# should probably have a slicker way of handling this
    otp_string = otp_string[1:]

    return otp_string

def stateframe_otp(lines, filename=None):
    '''Writes out the stateframe list '''
    if filename is None:
        fname = 'stateframe_x.h'
    else: 
        fname = filename

    f = open(fname, 'w')

    #prepend a date and time, first get suffix, could be 'h', 'c', 'f90'
    if fname.endswith('.f90'):
        comment_char_a = '! '
        comment_char_b = ' !'
    elif fname.endswith('.h'):
        comment_char_a = '/*'
        comment_char_b = '*/'
    elif fname.endswith('.c'):
        comment_char_a = '/*'
        comment_char_b = '*/'
    else:
        comment_char_a = '# '
        comment_char_b = ' #'

    #Add comment with filename, date and time
    localtime = time.asctime(time.localtime(time.time()))
    line0 = comment_char_a+' FILE: '+fname+' Created: '+localtime+comment_char_b

    s = str(line0)+'\n'
    f.write(s)
    for ppp in lines:
        s = str(ppp)+'\n'
        f.write(s)

    f.close()
    
def c2f90_struct(lines):
    '''Takes an unnested  C structure definintion for the stateframe, and
    writes the appropriate f90 type declarations'''

    nlines = len(lines)
    # first strip all semicolons
    arrx1 = range(nlines)
    for j in range(nlines):
        abc = lines[j]
        fff = len(abc)
        if abc[fff-1] == ';':
            arrx1[j] = abc[:fff-1]
        else:
            arrx1[j] = abc

    # create a 2-d list with each element
    arrx2 = range(nlines)
    for j in range(nlines):
        arrx2[j] = str.split(arrx1[j])

    otp_string = [-1]
    saved_structname = 'oops'
    for j in range(nlines):
        if arrx2[j][0] == 'struct':
            if arrx2[j][2] == '{':
            # this is a structure definition
                tstring = 'type '+arrx2[j][1]
                otp_string = otp_string+[tstring]+['   sequence']
                saved_structname = arrx2[j][1]
            else:
            # this is a variable definition, there may be dimensions
                vname = c2f90_reverse_indices(arrx2[j][2])
                tstring = '   type ('+arrx2[j][1]+') '+vname
                otp_string = otp_string+[tstring]
        elif arrx2[j][0] == '}':
            #this ends the structure, we've saved the structure name
            tstring = 'end type '+saved_structname
            otp_string = otp_string+[tstring]
            saved_structname = 'oops'
        else:
            #just plain variables
            if arrx2[j][0] == 'unsigned':
                if arrx2[j][1] == 'char':
                    vtype = 'integer *1'
                    vname = arrx2[j][2]
                elif arrx2[j][1] == 'short':
                    vtype = 'integer *2'
                    vname = arrx2[j][3]
                elif arrx2[j][1] == 'long':
                    vtype = 'integer *8'
                    vname = arrx2[j][3]
                else:
                    vtype = 'integer *4'
                    vname = arrx2[j][2]
            elif arrx2[j][0] == 'short':
                vtype = 'integer *2'
                vname = arrx2[j][1]
            elif arrx2[j][0] == 'int':
                vtype = 'integer *4'
                vname = arrx2[j][1]
            elif arrx2[j][0] == 'long':
                vtype = 'integer *8'
                vname = arrx2[j][1]
            elif arrx2[j][0] == 'float':
                vtype = 'real'
                vname = arrx2[j][1]
            elif arrx2[j][0] == 'double':
                vtype = 'real *8'
                vname = arrx2[j][1]
            elif arrx2[j][0] == 'char':
                vname = arrx2[j][1]
                #extract the array size (there should always be one)
                x0 = str.find(vname, '[')
                y0 = str.find(vname, ']')
                arrsz = vname[x0+1:y0]
                vname = vname[:x0]
                vtype = 'character'+'('+arrsz+')'
            #deal with any array dimensions, for non char vars
            vname = c2f90_reverse_indices(vname)
            tstring = '   '+vtype+' '+vname
            otp_string = otp_string+[tstring]
    #Should cover everything
    otp_string = otp_string[1:]
    return otp_string


def c_reverse_indices(a):
    '''Takes a c array, e.g., x[1][2][3], reverses indices, e.g., x[1][2][3]'''

    # find brackets
    c = str.count(a, '[')
    if c > 0:
        x0 = str.find(a, '[')
        y0 = str.find(a, ']')
        x1 = [x0]
        y1 = [y0]
        for j in range(c-1):
            x0 = str.find(a[x0+1:],'[')+x0+1
            x1 = x1+[x0]
            y0 = str.find(a[y0+1:],']')+y0+1
            y1 = y1+[y0]

        #reverse indices
        oxpt = a[:x1[0]]+'['
        for j in range(c):
            jj = c-j-1
            djj = a[x1[jj]+1:y1[jj]]
            oxpt = oxpt+djj
            if j == (c-1):
                oxpt = oxpt+']'
            else:
                oxpt = oxpt+']['
            
    else:
        oxpt = a

    return oxpt

def c2f90_reverse_indices(a):
    '''Takes a c array, e.g., x[1][2][3], returns a fortran array,
    e.g., x(3,2,1)'''

    # find brackets
    c = str.count(a, '[')
    if c > 0:
        x0 = str.find(a, '[')
        y0 = str.find(a, ']')
        x1 = [x0]
        y1 = [y0]
        for j in range(c-1):
            x0 = str.find(a[x0+1:],'[')+x0+1
            x1 = x1+[x0]
            y0 = str.find(a[y0+1:],']')+y0+1
            y1 = y1+[y0]

        #reverse indices, use parentheses and commas
        oxpt = a[:x1[0]]+'('
        for j in range(c):
            jj = c-j-1
            djj = a[x1[jj]+1:y1[jj]]
            oxpt = oxpt+djj
            if j == (c-1):
                oxpt = oxpt+')'
            else:
                oxpt = oxpt+','
            
    else:
        oxpt = a

    return oxpt

def c_append_stateframecomm(lines):
    '''Appends the definition of the stateframe global external structure for the stateframe'''

    nlines = len(lines)
    lines_out = range(nlines)
    for j in range(nlines):
        lines_out[j] = lines[j]

    # replace the last line here with a line that doesn't include the variable name
    lines_out[nlines-1] = '   };'
    lines_out = lines_out+['extern struct xstateframecomm {', '   struct xstateframe stateframe;', '   } stateframecomm_;']

    return lines_out

def c2f90_append_stateframecomm(lines):
    '''Appends the definition of the stateframe f90 common block'''

    nlines = len(lines)
    lines_out = range(nlines)
    for j in range(nlines):
        lines_out[j] = lines[j]

    lines_out = lines_out+['type (xstateframe) stateframe', 'common/stateframecomm/stateframe']

    return lines_out

def c_add_packed_attribute(lines):
    '''Adds __attribute__ ((__packed__)) to C structure definitions'''
    
    ext_string = "struct __attribute__ ((__packed__))"
    nlines = len(lines)
    # create a 2-d list for each element
    arrx2 = range(nlines)
    arrx1 = range(nlines)
    for j in range(nlines):
        arrx1[j] = lines[j]
        arrx1[j]=arrx1[j].replace('struct x', 'struct y', 1)
        arrx2[j] = str.split(arrx1[j])

    # Now for each structure definition, identified by arrx2[j][0] = "struct" and
    # arrx2[j][2] = "{"    
    for j in range(nlines):
        if arrx2[j][0] == 'struct' and arrx2[j][2] == '{':
            arrx1[j]=arrx1[j].replace('struct', ext_string)
    # you need to renamethe structure type too, instead of beginning all struct definitions with 'x', use 'y'
    return arrx1
    
def c_struct_varnames(lines, var0):
    '''Takes a C structure definintion for the stateframe, and creates
    a list of variable names and types that can be used in C code'''

    # first strip semicolons
    nlines = len(lines)
    arrx1 = range(nlines)
    for j in range(nlines):
        abc = lines[j]
        fff = len(abc)
        if abc[fff-1] == ';':
            arrx1[j] = abc[:fff-1]
        else:
            arrx1[j] = abc

    # create 2-d lists with each element
    arrx0 = range(nlines)
    arrx2 = range(nlines)
    for j in range(nlines):
        arrx0[j] = str.split(arrx1[j], ' ')
        arrx2[j] = str.split(arrx1[j])

    # str.split with a ' ' leaves null strings for spaces, we can
    # count to get indentation
    # for each line, find the first value and 
    # create a copy of arrx0 without the spaces in front
    level = range(nlines)
    first_x = range(nlines)
    for j in range(nlines):
        a = 0
        while arrx0[j][a] == '':
            a = a+1
        level[j] = a/3
        first_x[j] = arrx0[j][a]

    # find the start of the structures
    xxx = [-1]
    for j in range(nlines):
        if first_x[j] == 'struct':
            xxx = xxx+[j]

    xxx = xxx[1:]
        
    nstruct = len(xxx)
    xxx1 = range(nstruct)
    for j in range(nstruct):
        xxx1[j] = nlines-1

    nest_level = range(nstruct)
    strname = range(nstruct)
    # Now get the name of the structures
    for j in range(nstruct):
        j1 = xxx[j]
        nest_level[j] = level[j1]
        # here i need to find the next point where first_x is '}' and level is nest_level[j]+1
        while j1 < (nlines-1) and not(first_x[j1] == '}' and level[j1] == nest_level[j]+1):
            j1 = j1+1
        xxx1[j] = j1
        strname[j] = arrx2[xxx1[j]][1]

    # Insert var0 variable name into the strname[0] spot
    strname[0] = var0
    # Here create a list which for each line has either zero,
    # or for the positions xxx of the structure starts,
    # the value of the appropriate xxx1, the structure end
    stest = range(nlines)
    for j in range(nlines):
        stest[j] = 0

    for j in range(nstruct):
       for i in range(nlines):
           if i == xxx[j]:
               stest[i] = xxx1[j]

    # step through the structures, and just write the variable names
    maxnl = max(nest_level)
    arrx10 = range(nlines)
    for j in range(nlines):
        arrx10[j] = range(maxnl+1)
        for k in range(maxnl+1):
            arrx10[j][k] = ''

    for j in range(nstruct):
        # can't do: arrx10[xxx[j]:xxx1[j]][nest_level[j]] = strname[j]
        nxxx1 = xxx1[j]-xxx[j]+1
        for i in range(nxxx1):
            arrx10[xxx[j]+i][nest_level[j]] = strname[j]

    # contract the arrx10 array
    sn1plus = range(nlines)
    for j in range(nlines):
        sn1plus[j] = ''
        if level[j] > 0:
            for k in range(level[j]):
                sn1plus[j] = sn1plus[j]+arrx10[j][k]+'.'
        else:
            sn1plus[j] = arrx10[j][0]+'.'

    otpp = [-1]
    otpp2 = [-1]
    for j in range(nlines):
        if first_x[j] != 'struct' and first_x[j] != '}':
            xxlen = len(arrx2[j])-1
            otpp = otpp+[sn1plus[j]+arrx2[j][xxlen]]
            typ = ''
            for k in range(xxlen):
                typ = typ+arrx2[j][k]
            otpp2 = otpp2+[typ]

    otpp = otpp[1:]
    otpp2 = otpp2[1:]

    return otpp, otpp2

def c_vareqvar(a, b):
    '''given a variable name, set it equal to itself in C syntax. Vars 
    a and b are assumed to be the same structure with different names'''
    #find dimensionality if any
    # find brackets
    loopvars = ['i', 'j', 'k', 'l', 'm', 'n']
    c = str.count(a, '[')
    if c > 0:
        x0 = str.find(a, '[')
        y0 = str.find(a, ']')
        x1 = [x0]
        y1 = [y0]
        for j in range(c-1):
            x0 = str.find(a[x0+1:],'[')+x0+1
            x1 = x1+[x0]
            y0 = str.find(a[y0+1:],']')+y0+1
            y1 = y1+[y0]

        # find dimensions
        djj = range(c)
        for j in range(c):
            djj[j] = a[x1[j]+1:y1[j]]

        otp = [-1]
        # set up loops
        for j in range(c):
            test_str = 'for('+loopvars[j]+'=0;'+loopvars[j]+'<'+djj[j]+';'+loopvars[j]+'++)'
            otp = otp+[test_str, '{']

        # set variable = variable
        a1 = a
        b1 = b
        for j in range(c):
            djj_brackets = '['+djj[j]+']'
            loopvarsj_brackets = '['+loopvars[j]+']'
            a1 = a1.replace(djj_brackets, loopvarsj_brackets, 1)
            b1 = b1.replace(djj_brackets, loopvarsj_brackets, 1)

        # close loops
        otp = otp+[a1+' = '+b1+';']
        for j in range(c):
            otp = otp+['}']
        otp = otp[1:]

    else:
        otp = [a+' = '+b+';']

    return otp

def c_var_write(a, typex, funit):
    '''given a variable name and type, genereate a c fprintf statement for output, input funit is the file descriptor '''

    print typex
    print typex == 'double'

    #Need a format code
    if typex == 'int':
        fmt = '%i'+'/'+'n'
    elif typex == 'shortint':
        fmt = '%i'+'/'+'n'
    elif typex == 'unsignedint':
        fmt = '%i'+'/'+'n'
    elif typex == 'unsignedshortint':
        fmt = '%i'+'/'+'n'
    elif typex == 'longint':
        fmt = '%i'+'/'+'n'
    elif typex == 'unsignedlongint':
        fmt = '%i'+'/'+'n'
    elif typex == 'double':
        fmt = '%f'+'/'+'n'
    elif typex == 'float':
        fmt = '%f'+'/'+'n'
    elif typex == 'char':
        fmt = '%c'+'/'+'n'
    elif typex == 'unsignedchar':
        fmt = '%i'+'/'+'n'
    else:
        fmt = '%i'+'/'+'n'

    #find dimensionality if any
    # find brackets

    loopvars = ['i', 'j', 'k', 'l', 'm', 'n']
    c = str.count(a, '[')
    if c > 0:
        x0 = str.find(a, '[')
        y0 = str.find(a, ']')
        x1 = [x0]
        y1 = [y0]
        for j in range(c-1):
            x0 = str.find(a[x0+1:],'[')+x0+1
            x1 = x1+[x0]
            y0 = str.find(a[y0+1:],']')+y0+1
            y1 = y1+[y0]

        # find dimensions
        djj = range(c)
        for j in range(c):
            djj[j] = a[x1[j]+1:y1[j]]

        otp = [-1]
        # set up loops
        for j in range(c):
            test_str = 'for('+loopvars[j]+'=0;'+loopvars[j]+'<'+djj[j]+';'+loopvars[j]+'++)'
            otp = otp+[test_str, '{']

        # Make fprintf statements

        a1 = a
        for j in range(c):
            djj_brackets = '['+djj[j]+']'
            loopvarsj_brackets = '['+loopvars[j]+']'
            a1 = a1.replace(djj_brackets, loopvarsj_brackets, 1)

        fprintf_string = 'fprintf('+funit+','+'"'+fmt+'", '+a1+');'

        # close loops
        otp = otp+[fprintf_string]
        for j in range(c):
            otp = otp+['}']
        otp = otp[1:]

    else:
        otp = 'fprintf('+funit+','+'"'+fmt+'", '+a+');'

    return otp

def c2f90_var_write(a, typex, funit):
    '''given a variable name and type, genereate a fortran write statement for output, input funit is the file unit number '''

    #find dimensionality if any
    # find brackets

    loopvars = ['i', 'j', 'k', 'l', 'm', 'n']
    c = str.count(a, '[')
    if c > 0:
        x0 = str.find(a, '[')
        y0 = str.find(a, ']')
        x1 = [x0]
        y1 = [y0]
        for j in range(c-1):
            x0 = str.find(a[x0+1:],'[')+x0+1
            x1 = x1+[x0]
            y0 = str.find(a[y0+1:],']')+y0+1
            y1 = y1+[y0]

        # find dimensions
        djj = range(c)
        for j in range(c):
            djj[j] = a[x1[j]+1:y1[j]]

        otp = [-1]
        # set up loops
        for j in range(c):
            test_str = 'do '+loopvars[j]+'= 1, '+djj[j]
            otp = otp+[test_str, '{']

        # Make fprintf statements

        a1 = a
        for j in range(c):
            djj_brackets = '['+djj[j]+']'
            loopvarsj_brackets = '['+loopvars[j]+']'
            a1 = a1.replace(djj_brackets, loopvarsj_brackets, 1)

        fprintf_string = 'fprintf('+funit+','+'"'+fmt+'", '+a1+');'

        # close loops
        otp = otp+[fprintf_string]
        for j in range(c):
            otp = otp+['}']
        otp = otp[1:]

    else:
        otp = 'fprintf('+funit+','+'"'+fmt+'", '+a+');'

    return otp

def stateframe_eq_stateframe(lines, var1, var2):
    '''Writes a C subroutine that takes the input nested C structure (output of
    c_xml_read_stateframe), and sets the padded stateframe structure var1
    equal to the unpadded one. var1 and 2 are assumed to be global stateframe
    variables and are not defined in this subroutine'''

    otp = ['#include        "Stateframe_packed.h"', '#include        "DPPstateframe_ext.h"']
    otp = otp+['/* Function Prototype */', 'void stateframe_eq_stateframe(void);']
    otp = otp+['/* Function Definition */', 'void stateframe_eq_stateframe(void)', '{']
    otp = otp+['int i,j,k,l,m; /*loop variables*/']
    
    var1x, typ1x = c_struct_varnames(lines, var1)
    var2x, typ1x = c_struct_varnames(lines, var2)

    for j in range(len(var1x)):
        otp = otp+c_vareqvar(var1x[j], var2x[j])

    otp = otp+['}', '/* End of function */']

    return otp

def write_stateframe_files():
    '''Reads the stateframe XML file, and writes out the structure definitions to be used by the DPP'''

    # Read the xml definition, 
    lines0, version0 = c_xml_read_stateframe()
    stateframe_otp(lines0, filename='Stateframe0.h')
    
    # Unnest C structure, and write files with packed attribute and
    #  with global external definition
    lines1 = c_struct_unnest(lines0)
    lines11 = c_add_packed_attribute(lines1)
    #Add the version number to the packed stateframe definition, where I hope it;ll get picked up everywhere
    lines11a = ["#define SAVED_STATEFRAME_VERSION "+version0+" /* Stateframe version number */"]+lines11
    stateframe_otp(lines11a, filename='Stateframe_packed.h')
    lines2 = c_append_stateframecomm(lines1)
    stateframe_otp(lines2, filename='DPPstateframe_ext.h')

    # Write FORTRAN type definition
    flines = c2f90_struct(lines1)
    flines1 = c2f90_append_stateframecomm(flines)
    stateframe_otp(flines1, filename='DPPstateframe.f90')

    # Write C code for the process that sets the padded structure 
    # equal to the packed structure.
    clines = stateframe_eq_stateframe(lines0, 'stateframecomm_.stateframe', 'Stateframe')
    stateframe_otp(clines, 'stateframe_eq_stateframe.c')

def write_scanheader_files():
    '''Same sort of treatment for the scan header'''
    sclines0, scversion0 = c_xml_read_scanheader()
    stateframe_otp(sclines0, filename='Scanheader0.h')

    # Unnest structure
    sclines1 = c_struct_unnest(sclines0)

    # Add packed attribute and version number
    sclines11 = c_add_packed_attribute(sclines1)
    sclines11a = ["#define SAVED_SCANHEADER_VERSION "+scversion0+" /* Scanheader version number */"]+sclines11
    stateframe_otp(sclines11a, filename = 'Scanheader_packed.h')

    # Add external structure definition
    sclines2 = c_append_scanheadercomm(sclines1)
    stateframe_otp(sclines2, filename = 'DPPscanheader_ext.h')

    # DO fortran output
    fsclines = c2f90_struct(sclines1)
    fsclines1 = c2f90_append_scanheadercomm(fsclines)
    stateframe_otp(fsclines1, filename = 'DPPscanheader.f90')

    # write C program to set packed to padded, tag by tag
    csclines = scanheader_eq_scanheader(sclines0, 'scanheadercomm_.scanheader', 'Scan_Header')
    stateframe_otp(csclines, 'scanheader_eq_scanheader.c')


def c_read_scanheader(filename=None):
    '''Reads the scanheader structure definition from the file "Scanheader_packed.h"'''

    if filename is None:
        f = open('Scanheader_packed.h', 'r')
    else: 
        f = open(filename)

    lines = f.readlines()

    #remove the \n
    nlines = len(lines)
    for j in range(nlines):
        temp_l = lines[j]
        len_l = len(temp_l)
        temp_l = temp_l[:len_l-2]
        lines[j] = temp_l

    f.close()

    return lines

def c_append_scanheadercomm(lines):
    '''Takes the packed scan header definition and returns a list for
    the padded scan header structure definition in C'''

    nlines = len(lines)
    lines_out = range(nlines)
    for j in range(nlines):
        lines_out[j] = lines[j]

    lines_out[nlines-1] = '   };'
    lines_out = lines_out+['extern struct xscanheadercomm {', 'struct xscan_header scanheader;', '} scanheadercomm_;']

    return lines_out

def c2f90_append_scanheadercomm(lines):
    '''Takes the packed scan header definition and returns a list for
    the padded scan header structure definition in FORTRAN'''



    nlines = len(lines)
    lines_out = range(nlines)
    for j in range(nlines):
        lines_out[j] = lines[j]

    lines_out = lines_out+['type (xscan_header) scanheader', 'common/scanheadercomm/scanheader']
    return lines_out

def scanheader_eq_scanheader(lines, var1, var2):
    '''Writes a C subroutine that takes the input nested C structure (output of
    c_xml_read_scanheader), and sets the padded scanheader structure var1
    equal to the unpadded one. var1 and 2 are assumed to be global scanheader
    variables and are not defined in this subroutine'''

    otp = ['#include        "Scanheader_packed.h"', '#include        "DPPscanheader_ext.h"']
    otp = otp+['/* Function Prototype */', 'void scanheader_eq_scanheader(void);']
    otp = otp+['/* Function Definition */', 'void scanheader_eq_scanheader(void)', '{']
    otp = otp+['int i,j,k,l,m; /*loop variables*/']
    
    var1x, typ1x = c_struct_varnames(lines, var1)
    var2x, typ1x = c_struct_varnames(lines, var2)

    for j in range(len(var1x)):
        otp = otp+c_vareqvar(var1x[j], var2x[j])

    otp = otp+['}', '/* End of function */']

    return otp


