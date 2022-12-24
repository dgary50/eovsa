# Routines for finding trajectory information (after using *.trj file to search)
# These were written when Antenna 3 had a very large offset (around 4 degrees), and
# a .trj file was used to search a square pattern in the sky to help location the
# Sun.
#
#  Written 2022-09-27  DG

def find_traj(trange):
    ''' Searches SQL for any TRAJ-FILE and TRAJ-ON commands within the given timerange,
        and prints some useful information to the terminal.  Specifically, it prints 
        the times of the TRAJ-FILE and TRAJ-ON commands and the current RA and Dec when 
        the TRAJ-ON command was issued.
    '''
    from util import Time
    import dbutil as db
    cursor = db.get_cursor()
    t1, t2 = trange.lv.astype(int)
    query = "select Timestamp,Sche_Task from fV66_vD1 where Timestamp between "+str(t1)+" and "+str(t2)
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        for i,task in enumerate(data['Sche_Task']):
            if task.find('TRAJ-FILE') != -1:
                print task, Time(data['Timestamp'][i],format='lv').iso
            elif task.find('TRAJ-ON') != -1:
                t1 = data['Timestamp'][i]
                print task, Time(t1,format='lv').iso
                sqldict = db.get_dbrecs(cursor,dimension=15,timestamp=t1,nrecs=1)
                print 'RA (degrees)',sqldict['Ante_Cont_RAVirtualAxis']/10000.
                print 'Dec (degrees)',sqldict['Ante_Cont_DecVirtualAxis']/10000.
    else:
        print 'Error in SQL query',msg
        
        
def trajoff2xelel(radeg, decdeg, t0, raoffuser, decoffuser):
    ''' Converts a measured RA and Dec offset (in "User" units of 1/10000th of a degree) 
        to XEL and EL offsets printed to the terminal, also in "User" units.
        
        Inputs are the RA and Dec, in degrees, the Time() object for the time of the
        measurement, and the RA and Dec offsets in "User" units.
    '''
    from solpnt import dradec2dazel
    import numpy as np
    dtor = np.pi/180.
    user2rad = dtor/10000.
    ra = radeg*dtor
    dec = decdeg*dtor
    cosdec = np.cos(dec)
    xel, el = dradec2dazel(radeg*dtor,decdeg*dtor, t0, raoffuser*user2rad/cosdec,decoffuser*user2rad)
    print 'Sun Offset in User units:', int(xel/user2rad), int(el/user2rad)