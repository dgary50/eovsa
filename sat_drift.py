#
#   Routine for creating a "drift" track table.
# 
# History:
#   2015-Oct-27  DG
#      First written

import urllib2
from util import Time
from ftplib import FTP
import time

def sat_drift(t,dir='RA',dist=1,dt=60.):
    ''' Doctor a geosat_tab.radec file to cause an antenna
        to move dist degrees in dt seconds, in direction 'RA' or 'Dec'
        Defaults are 'RA' at 1 degree/minute.

        Usage: sat_drift('2015-10-26 22:00:00'[,dir='Dec'][,dist=2][,dt=180.])

        The above options will cause 2-degree drift in 3 minutes in Dec
    '''
    userpass = 'admin:observer@'
    f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/geosat_tab.radec',timeout=0.5)
    lines = f.readlines()
    f.close()
    mjd2 = Time(t).mjd # Middle time is given time
    for i,line in enumerate(lines):
        ra, da, mjd, ta = line.split()
        mjda = int(mjd) + float(ta)/86400000.
        if mjda > mjd2:
            rb, db, mjd, tb = lines[i-1].split()
            mjdb = int(mjd) + float(tb)/86400000.
            print i
            break

    dtim = dt/86400.  # Convert to days

    mjd0 = mjd2 - dtim
    mjd1 = mjd2 - dtim/2
    mjd3 = mjd2 + dtim/2
    mjd4 = mjd2 + dtim

    if dir == 'RA':
        dra = dist*10000.
        ddec = 0.
    elif dir == 'Dec':
        dra = 0.
        ddec = dist*10000.

    ra = float(ra)
    rb = float(rb)
    ta = float(ta)
    tb = float(tb)
    da = float(da)
    db = float(db)

    rb1= (ra - rb)*(mjd0 - mjdb)/(mjda - mjdb) + rb  # Time 0 with no position offset
    r0 = (ra - rb)*(mjd0 - mjdb)/(mjda - mjdb) + rb - dra
    r1 = (ra - rb)*(mjd1 - mjdb)/(mjda - mjdb) + rb - dra/2
    r2 = (ra - rb)*(mjd2 - mjdb)/(mjda - mjdb) + rb
    r3 = (ra - rb)*(mjd3 - mjdb)/(mjda - mjdb) + rb + dra/2
    r4 = (ra - rb)*(mjd4 - mjdb)/(mjda - mjdb) + rb + dra
    ra1= (ra - rb)*(mjd4 - mjdb)/(mjda - mjdb) + rb  # Time 4 with no position offset

    db1= (da - db)*(mjd0 - mjdb)/(mjda - mjdb) + db # Time 0 with no position offset
    d0 = (da - db)*(mjd0 - mjdb)/(mjda - mjdb) + db - ddec
    d1 = (da - db)*(mjd1 - mjdb)/(mjda - mjdb) + db - ddec/2.
    d2 = (da - db)*(mjd2 - mjdb)/(mjda - mjdb) + db
    d3 = (da - db)*(mjd3 - mjdb)/(mjda - mjdb) + db + ddec/2.
    d4 = (da - db)*(mjd4 - mjdb)/(mjda - mjdb) + db + ddec
    da1= (da - db)*(mjd4 - mjdb)/(mjda - mjdb) + db # Time 4 with no position offset

    tbl = ''
    for line in lines[:i]:
        tbl += line.rstrip()+'\n'
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(rb1), int(d2), int(mjd0), int((mjd0 % 1)*86400000.+ 0.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(rb1), int(d2), int(mjd0), int((mjd0 % 1)*86400000.+10.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(r0),  int(d0), int(mjd0), int((mjd0 % 1)*86400000.+20.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(r1),  int(d1), int(mjd1), int((mjd1 % 1)*86400000.+ 0.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(r2),  int(d2), int(mjd2), int((mjd2 % 1)*86400000.+ 0.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(r3),  int(d3), int(mjd3), int((mjd3 % 1)*86400000.+ 0.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(r4),  int(d4), int(mjd4), int((mjd4 % 1)*86400000.-19.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(ra1), int(da), int(mjda), int((mjda % 1)*86400000.- 9.5))
    tbl += '{:>7d} {:>7d} {:5} {:>8d}\n'.format(int(ra1), int(da), int(mjda), int((mjda % 1)*86400000.+ 0.5))
    for line in lines[i:]:
       tbl += line.rstrip()+'\n'

    print tbl
    fname = 'geodrift_tab.radec'
    f = open('/tmp/'+fname,'w')
    f.write(tbl)
    f.close()
    time.sleep(0.01)
    # Send tracktable file to acc
    f = open('/tmp/'+fname,'r')
    acc = FTP('acc.solar.pvt')
    acc.login('admin','observer')
    acc.cwd('parm')
    acc.storlines('STOR '+fname,f)
    acc.close()
    f.close()



