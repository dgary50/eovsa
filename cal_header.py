#
# Routines for implementing XML descriptions of calibration data.
# 
# History:
#   2015-Mar-29  DG
#     First written (for Solpnt).
#   2015-Mar-30  DG
#     Added FGHz and Poln arrays to Solpnt
#   2015-Mar-31  DG
#     Added DCM_base_attn().  Changed name Solpnt to TPcal, and changed 
#     routine names to end in "...2xml"  Also added cal_types dictionary.
#   2015-Apr-01  DG
#     Forgot to commit changes to database!  Added commit at bottom of
#     send_xml2sql().
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#
import struct, util

def cal_types():
    ''' Routine that defines all of the basic "long-term" calibration
        /information types as a dictionary.  These types and descriptions 
        will be written into the Description field of the aBin table.
        A new type can be added at the end--there is no significance to
        the type number--it is just a unique ordinal.
    '''
    return {'1':['Total power calibration (output of SOLPNTCAL)','TPcal2xml'],
            '2':['DCM master base attenuation table [units=dB]','DCM_master_base_attn2xml'],
            '3':['DCM base attenuation table [units=dB]','DCM_base_attn2xml'],
            '4':['Delay centers [units=ns]','dlacen2xml']}

def str2bin(string):
    return struct.pack(str(len(string)+1)+'s',string+'\n')

def TPcal2xml():
    ''' Writes the XML description of the total power calibration binary
        data (SOLPNTCAL result).  Returns a binary representation of the
        text file, for putting into the SQL database.  The version number
        must be incremented each time there is a change to the structure 
        of this header.
    '''
    version = 1.0
    
    buf = ''
    buf += str2bin('<Cluster>')
    buf += str2bin('<Name>TPcal</Name>')
    buf += str2bin('<NumElts>5</NumElts>')
    
    # Timestamp (double) [s, in LabVIEW format]
    # Start time of SOLPNT observation on which calibration is based
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Timestamp</Name>')
    buf += str2bin('<Val></Val>')
    buf += str2bin('</DBL>')
    
    # Version of this XML file.  This number should be incremented each
    # time there is a change to the structure of this file.
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Version</Name>')
    buf += str2bin('<Val>'+str(version)+'</Val>')
    buf += str2bin('</DBL>')

    # Array of frequencies in GHz (448)
    buf += str2bin('<Array>')
    buf += str2bin('<Name>FGHz</Name>')
    buf += str2bin('<Dimsize>448</Dimsize>\n<SGL>\n<Name></Name>\n<Val></Val>\n</SGL>')
    buf += str2bin('</Array>')

    # Array of Poln (2)
    # Polarization list (Miriad definition) (signed int)
    #     1: Stokes I
    #     2: Stokes Q
    #     3: Stokes U
    #     4: Stokes V
    #    -1: Circular RR
    #    -2: Circular LL
    #    -3: Circular RL
    #    -4: Circular LR
    #    -5: Linear XX
    #    -6: Linear YY
    #    -7: Linear XY
    #    -8: Linear YX
    #     0: Not used
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Poln</Name>')
    buf += str2bin('<Dimsize>2</Dimsize>\n<I32>\n<Name></Name>\n<Val></Val>\n</I32>')
    buf += str2bin('</Array>')

    # Array of clusters for each antenna (13 since for 2.1-m ants only)
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Antenna</Name>')
    buf += str2bin('<Dimsize>13</Dimsize>')

    # Cluster containing information for one antenna
    buf += str2bin('<Cluster>')
    buf += str2bin('<Name></Name>')
    buf += str2bin('<NumElts>2</NumElts>')

    # Calibration factors (448 x 2) = nfreq x npol
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Calfac</Name>')
    buf += str2bin('<Dimsize>448</Dimsize><Dimsize>2</Dimsize>\n<SGL>\n<Name></Name>\n<Val></Val>\n</SGL>')
    buf += str2bin('</Array>')

    # Offsun values (448 x 2) = nfreq x npol
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Offsun</Name>')
    buf += str2bin('<Dimsize>448</Dimsize><Dimsize>2</Dimsize>\n<SGL>\n<Name></Name>\n<Val></Val>\n</SGL>')
    buf += str2bin('</Array>')

    # End cluster
    buf += str2bin('</Cluster>')   # End Calinfo cluster
    buf += str2bin('</Array>')     # End Antenna array
    buf += str2bin('</Cluster>')   # End TPcal cluster

    return buf

def DCM_master_base_attn2xml():
    ''' Writes the XML description of the DCM master base attenuation 
        table (created by pcapture.py).  Returns a binary representation 
        of the text file, for putting into the SQL database.  The version 
        number must be incremented each time there is a change to the 
        structure of this header.
    '''
    version = 1.0
    
    buf = ''
    buf += str2bin('<Cluster>')
    buf += str2bin('<Name>DCMMasterBaseAttn</Name>')
    buf += str2bin('<NumElts>5</NumElts>')
    
    # Timestamp (double) [s, in LabVIEW format]
    # Time of creation of the table (close to packet capture time)
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Timestamp</Name>')
    buf += str2bin('<Val></Val>')
    buf += str2bin('</DBL>')
    
    # Version of this XML file.  This number should be incremented each
    # time there is a change to the structure of this file.
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Version</Name>')
    buf += str2bin('<Val>'+str(version)+'</Val>')
    buf += str2bin('</DBL>')

    # List of antennas (15), with antenna number (1-15) if used, 0 if not.
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Antennas</Name>')
    buf += str2bin('<Dimsize>15</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # List of bands (34), with band number (1-34) if used, 0 if not.
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Bands</Name>')
    buf += str2bin('<Dimsize>34</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # Array of base attenuations [dB] (34 x 30).  Attenuations for unmeasured
    # antennas and/or bands are set to nominal value of 10 dB.  Values are
    # ordered as Ant1x, Ant1y, Ant2x, Ant2y, ..., Ant15x, Ant15y
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Attenuation</Name>')
    buf += str2bin('<Dimsize>30</Dimsize><Dimsize>34</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # End cluster
    buf += str2bin('</Cluster>')   # End DCMMasterBaseAttn cluster

    return buf

def DCM_base_attn2xml():
    ''' Writes the XML description of the DCM base attenuation table 
        (derived from the DCM master base attenuation table and the
        current frequency sequence.  Returns a binary representation of the
        text file, for putting into the SQL database.  The version number
        must be incremented each time there is a change to the structure 
        of this header.
    '''
    version = 1.0
    
    buf = ''
    buf += str2bin('<Cluster>')
    buf += str2bin('<Name>DCMBaseAttn</Name>')
    buf += str2bin('<NumElts>5</NumElts>')
    
    # Timestamp (double) [s, in LabVIEW format]
    # Time of creation of the table (close to packet capture time)
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Timestamp</Name>')
    buf += str2bin('<Val></Val>')
    buf += str2bin('</DBL>')
    
    # Version of this XML file.  This number should be incremented each
    # time there is a change to the structure of this file.
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Version</Name>')
    buf += str2bin('<Val>'+str(version)+'</Val>')
    buf += str2bin('</DBL>')

    # List of antennas (15), with antenna number (1-15) if used, 0 if not.
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Antennas</Name>')
    buf += str2bin('<Dimsize>15</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # List of bands in 50-slot frequency sequence (50), with band number (1-34).
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Bands</Name>')
    buf += str2bin('<Dimsize>50</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # Array of base attenuations [dB] (50 x 30).  Attenuations for unmeasured
    # antennas and/or bands are set to nominal value of 10 dB.  Values are
    # ordered as Ant1x, Ant1y, Ant2x, Ant2y, ..., Ant15x, Ant15y
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Attenuation</Name>')
    buf += str2bin('<Dimsize>30</Dimsize><Dimsize>50</Dimsize>\n<U16>\n<Name></Name>\n<Val></Val>\n</U16>')
    buf += str2bin('</Array>')

    # End cluster
    buf += str2bin('</Cluster>')   # End DCMBaseAttn cluster

    return buf

def dlacen2xml():
    ''' Writes the XML description of the Delay Centers table (currently
        created by hand).  Returns a binary representation of the xml
        text file, for putting into the SQL database.  The version number
        must be incremented each time there is a change to the structure 
        of this header.
    '''
    version = 1.0
    
    buf = ''
    buf += str2bin('<Cluster>')
    buf += str2bin('<Name>Delaycenters</Name>')
    buf += str2bin('<NumElts>3</NumElts>')
    
    # Timestamp (double) [s, in LabVIEW format]
    # Time of creation of the table (precise time not critical)
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Timestamp</Name>')
    buf += str2bin('<Val></Val>')
    buf += str2bin('</DBL>')
    
    # Version of this XML file.  This number should be incremented each
    # time there is a change to the structure of this file.
    buf += str2bin('<DBL>')
    buf += str2bin('<Name>Version</Name>')
    buf += str2bin('<Val>'+str(version)+'</Val>')
    buf += str2bin('</DBL>')

    # List of delay centers [ns] (16).
    buf += str2bin('<Array>')
    buf += str2bin('<Name>Delaycen_ns</Name>')
    buf += str2bin('<Dimsize>16</Dimsize>\n<SGL>\n<Name></Name>\n<Val></Val>\n</SGL>')
    buf += str2bin('</Array>')

    # End cluster
    buf += str2bin('</Cluster>')   # End Delaycenters cluster

    return buf

def send_xml2sql():
    ''' Routine to send any changed calibration xml definitions to the 
        SQL Server.  The latest definition (if any) for a given type is
        checked to see if the version matches.  If not, an update is 
        stored.  This routine will typically be run whenever a definition
        is added or changed.
    '''
    import dbutil, read_xml2, sys
    t = util.Time.now()
    timestamp = t.lv  # Current time as LabVIEW timestamp
    cursor = dbutil.get_cursor()
    typdict = cal_types()
    for key in typdict.keys():
        #print 'Working on',typdict[key][0]
        # Execute the code to create the xml description for this key
        exec 'buf = '+typdict[key][1]+'()'
        # Resulting buf must be written to a temporary file and reread
        # by xml_ptrs()
        f = open('/tmp/tmp.xml','wb')
        f.write(buf)
        f.close()
        mydict, xmlver = read_xml2.xml_ptrs('/tmp/tmp.xml')
        defn_version = float(key)+xmlver/10.  # Version number expected
        # Retrieve most recent key.0 record and check its version against the expected one
        query = 'select top 1 * from abin where Version = '+key+'.0 order by Timestamp desc'
        #print 'Executing query'
        outdict, msg = dbutil.do_query(cursor,query)
        #print msg
        if msg == 'Success':
            if len(outdict) == 0:
                # This type of xml file does not yet exist in the database, so mark it for adding
                add = True
            else:
                # There is one, so see if it agrees with the version
                buf = outdict['Bin'][0]   # Binary representation of xml file
                f = open('/tmp/tmp.xml','wb')
                f.write(buf)
                f.close()
                mydict, thisver = read_xml2.xml_ptrs('/tmp/tmp.xml')
                #print 'Versions:',float(key)+thisver/10.,defn_version
                if (float(key)+thisver/10.) == defn_version:
                    # This description is already there so we will skip it
                    add = False
                else:
                    add = True
        else:
            # Some kind of error occurred, so print the message and skip adding the xml description
            #print 'Query',query,'resulted in error message:',msg
            add = False
        if add:
            # This is either new or updated, so add the xml description
            # to the database
            #print 'Trying to add',typdict[key][0]
            try:
                cursor.execute('insert into aBin (Timestamp,Version,Description,Bin) values (?,?,?,?)',
                   timestamp, float(key), typdict[key][0], dbutil.stateframedef.pyodbc.Binary(buf))
                print typdict[key][0],'successfully added/updated to version',defn_version
            except:
                print 'Unknown error occurred in adding',typdict[key][0]
                print sys.exc_info()[1]
        else:
            print typdict[key][0],'version',defn_version,'already exists--not updated'
    cursor.commit()
    cursor.close()
            
            
            
            
