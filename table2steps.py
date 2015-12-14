import numpy as np

def table2steps(tbl,p=None):
    '''Converts an RADEC track table into a table in units of the number of motor microsteps 
       needed to point at the locations in the table.
    '''
    lat = 37.233170   # Latitude  of OVRO [deg]
    lng = -118.286953 # Longitude of OVRO [deg]
    dtor = np.pi/180. # Convert degrees to radians
    # Just a random set of pointing parameters, to exercise that code
    if p is None: p = [0]*9
    # Test parameters: p = [0.123, -0.42, -0.0013, 0.213, -0.001, -0.412, 0.44, -0.23, 0.0]
    xstepsiz = 15886. # Number of microsteps per degree in X
    ystepsiz = 15886. # Number of microsteps per degree in Y

    # Loop through lines in the table
    for line in tbl:
        # Extract values from line
        ra, dec, mjd, t = np.array(line.strip().split(),'int')
        # convert Dec to degrees
        dec = dec/10000.
        # Calculate days (float, includes time of day) since Jan. 1, 2000, at 12 UT
        day = mjd + t/86400000. - 51544.5
        # Convert RA at this time to HA [deg]
        ha = (18.697374558 + 24.06570982441908*day)*15 - ra/10000. + lng
        # Truncate HA to LHA between -180. and 180.
        lha = ha % 360.
        if lha > 180.: lha -= 360.

        # Convert angles to radians and calculate pointing offsets 
        # based on pointing parameters P
        rlat = lat*dtor
        rdec = dec*dtor
        rha = lha*dtor
        ho = p[0] - p[1]*np.cos(rlat)*np.sin(rha)/np.cos(rdec) + p[2]*np.tan(rdec) \
                  - p[3]/cos(rdec) + p[4]*sin(rha)*tan(rdec) - p[5]*cos(rha)*tan(rdec)
        do = p[4]*cos(rha) + p[5]*sin(rha) + p[6] - p[7]*(cos(rlat)*cos(rha)*sin(rdec) 
                  - sin(rlat)*cos(rdec))

        # Correct for offsets
        lha += ho
        dec += do

        # Convert to steps and print table
        x = int(lha*xstepsiz)
        y = int((dec - lat)*ystepsiz)
        print '{:8d} {:8d} {:5d} {:8d}'.format(x, y, mjd, t)

