from lxml import etree
import numpy as np
import copy
import struct
import urllib2

global indent
indent = 0

def c_handle_cluster(child):
    '''This element of the XML tree is the head of a Cluster.  Step through
       each element of the branch and return the keys, the empty dictionary, 
       and the fmt string.  This routine is reentrant'''
    global indent
    # Clusters have a name, a NumElts, and one or more objects
    c = list(child)
    if c[0].tag == "Name":
        sname = c[0].text  # Name of structure
        lines = []
        if sname:
            lines = ['struct x'+sname.lower()+' {'] # Start struct line list
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
            for dim in dims:
                line += str(dim)+','
            line = line[:len(line)-1]+'];'
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
            lines += [' '*3+'   } '+name+';']
            c.pop(0)
        else:
            ch = c[0]
            name, datatype = c_handle_item(ch)
            # Add a simple line entry with datatype followed by name of variable
            lines += [' '*3+datatype+' '+name+';']
            c.pop(0)
    return sname, lines

def c_handle_array(child):
    '''This element of the XML tree is an Array.  Step through the items of the
       array (which may contain clusters, other arrays, etc.) and return the keys,
       array, dimensions of the array, and fmt string.  This routine is reentrant.
    '''
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
    objecttype = c[0].tag
    dtype_dict = {'U8':'char','B8':'char','U16':'unsigned short int','U32':'unsigned int',
                  'I16':'short int','I32':'int','SGL':'float','DBL':'double'}
    datatype = dtype_dict.get(objecttype,'*')
    if objecttype == "Cluster":
        ch = list(c[0])
        name2, lines = c_handle_cluster(ch)
        if name2:
            name = name2
        datatype = 'struct x'+name.lower()+' {'
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

def c_xml_read(filename=None):
    '''Read an XML description of the stateframe and return lines for output
    to a C-language header file.  This header file can be included in a
    compilation to allow reading the binary data corresponding to the XML file
    and memory mapping it to the structure.
    '''
    if filename is None:
        f = urllib2.urlopen('ftp://acc.solar.pvt/ni-rt/startup/stateframe.xml')
    else: 
        f = open(filename)
    tree = etree.parse(f)
    f.close()

    root = tree.getroot()
    name, lines = c_handle_cluster(root)
    lines.append('   } '+name+';')
    return lines


