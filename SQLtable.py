#
# Routines for writing calibration data into SQL database.
# 
# History:
#   2015-Apr-02  DG
#     First written.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy
#

import struct, dbutil
from util import Time
import numpy as np

def TPcal(x, y, calfac, offsun):
    ''' Writes Total Power calibration factors and offsun IF level
        to SQL database (caltype = 1)
    '''
    # *******
    # Version has to be updated at same time as xml description in cal_header.py
    # *******
    version = 1.0  
    fghz = x['fghz']
    nf = len(fghz)
    tstamp = x['ut'][0] # Start time of SOLPNTCAL
    dims = calfac.shape
    buf = ''
    if nf == 448:
        buf = struct.pack('d',tstamp)
        buf += struct.pack('d',version)
        # Case of 448 frequencies only
        # Array dimension for frequency list
        buf += struct.pack('I',448)
        # Frequency list
        buf += struct.pack('448f',*fghz)
        # Polarization array (dimension, then two states--XX, YY)
        buf += struct.pack('Iii',*[2,-5,-6])
        # Array dimension for Antenna cluster (2.1 m ants only)
        buf += struct.pack('I',13)
        # Empty array for filling in for missing antennas
        empty = np.zeros(448,'float')
        for i in range(dims[2]):
            # Array dimensions for freq/poln for this antenna
            buf += struct.pack('2I',*[448,2])
            # Cal factors for the two polarizations
            buf += struct.pack('448f',*calfac[0,:,i])
            buf += struct.pack('448f',*calfac[1,:,i])
            # Array dimensions for freq/poln for this antenna
            buf += struct.pack('2I',*[448,2])
            # Offsun IF level for the two polarizations
            buf += struct.pack('448f',*offsun[0,:,i])
            buf += struct.pack('448f',*offsun[1,:,i])
        for i in range(dims[2],13):
            # Same as above for missing antennas
            buf += struct.pack('2I',*[448,2])
            buf += struct.pack('448f',*empty)
            buf += struct.pack('448f',*empty)
            buf += struct.pack('2I',*[448,2])
            buf += struct.pack('448f',*empty)
            buf += struct.pack('448f',*empty)
        t = Time.now()
        timestamp = t.lv
        cursor = dbutil.get_cursor()
        cursor.execute('insert into aBin (Timestamp,Version,Description,Bin) values (?,?,?,?)',
                   timestamp,1.0+version/10.,'Total Power Calibration',dbutil.stateframedef.pyodbc.Binary(buf))
        # *******
        #  NB! To retrieve these large binary data strings, one must declare text size on select, e.g.
        #         cursor.execute('set textsize 100000 select * from aBin where version = 1.1 ')
        #  where the given size is greater than the size desired.
        # *******
        
        # Temporarily store in disk file for checking format...
        #f = open('/tmp/tpcal.dat','wb')
        #f.write(buf)
        #f.close()

        cursor.commit()
        cursor.close()
