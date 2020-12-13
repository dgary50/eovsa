#
# EOVSA_MOUNTCAL
#  Routine to convert pointing measurements analyzed with calibration.py 
#  routines calpntanal() and calpnt_multi(), to the appropriate pointing
#  parameters used in the antennas.
#   
# History
#  2018-01-20  DG
#    Finally wrote a purely EOVSA version of mountcal that can handle all
#    of the 2.1-m antennas.  Needs work to add 27-m, due to broken nature
#    of 27-m pointing model.  I have been through this pretty carefully, and
#    it should be correct.
#  2018-01-30  DG
#    Added params2ants() function to send the new pointing parameters
#  2019-06-24  DG
#    Added Ant 14 broken pointing model (it was a very easy update to the
#    EQ mount type).
#  2020-11-29  DG
#    Added '.' for path in case of path == '' when creating plot.
#
from util import Time
import numpy as np
from eovsa_lst import eovsa_lst
dtor = np.pi/180.  # Converts degrees to radians
lat = 37.233170*np.pi/180
lng = -118.286953*np.pi/180

def rd_calpnt(filename):
    ''' Read and return contents of output of CALPNT or CALPNT2M observation
        Note that the "x" offsets are dRA and dAZ.  To apply these in subsequent
        routines, dRA is converted to dHA by inverting the sign, while dAZ is 
        converted to dXEL (angle on the sky) by multiplying by cos(EL).
    '''
    import dbutil
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    ants = lines[0].split('Ant ')[1:]
    antlist = []
    try:
        for ant in ants:
            antlist.append(int(ant[:2]))
    except:
        print 'Unrecognized format (header line) for',filename
        return None
    nants = len(ants)
    lines = lines[1:]
    nlines = len(lines)
    dra = np.zeros((nants,nlines),np.float)
    ddec = np.zeros((nants,nlines),np.float)
    ha = []
    dec = []
    source = []
    timstr = []
    for i,line in enumerate(lines):
        if len(line) < 9:
            dra = dra[:,:i]
            ddec = ddec[:,:i]
            break
        vals = line[9:].split()
        if len(vals) != nants*2 + 4:
            print 'Error reading line',i+2,'of',filename
            dra = dra[:,:i]
            ddec = ddec[:,:i]
            break
        else:
            try:
                source.append(line[:8])
                timstr.append(vals[0]+' '+vals[1])
                ha.append(vals[2])
                dec.append(vals[3])
                dra[:,i] = np.array(vals[4::2]).astype(float)
                ddec[:,i] = np.array(vals[5::2]).astype(float)
            except:
                print 'Error parsing line',i+2,'of',filename
    # Convert HA, Dec from degrees to radians, and then convert HA to RA
    ha = np.array(ha).astype(float)*np.pi/180.
    dec = np.array(dec).astype(float)*np.pi/180.
    times = Time(timstr)
    ra = np.zeros_like(ha)
    for i in range(len(ha)):
        ra[i] = eovsa_lst(times[i]) - ha[i]

    # Read pointing parameters from SQL database at time of first observation
    params_old = np.zeros((9,15),int)
    cursor = dbutil.get_cursor()
    timestamp = times[0].lv
    # Read stateframe data at time of first observation
    D15data = dbutil.get_dbrecs(cursor, dimension=15, timestamp=timestamp, nrecs=1)
    for p in range(9):
        params_old[p], = D15data['Ante_Cont_PointingCoefficient'+str(p+1)]
    params_old = params_old[:,np.array(antlist)-1]  # Pare down to only antennas in antlist

    return {'filename':filename, 'source':source, 'time':times, 'params_old':params_old, 'ra':ra, 
            'dec':dec, 'ha':ha, 'antlist':antlist, 'dra':dra, 'ddec':ddec}
            
def mntcal(indict):
    ''' Given an input dictionary of pointing coordinates for a single antenna,
        calculate pointing parameters according to antenna type.
        
        Input:
          indict   A dictionary with the following keys:
                     'ant'   Antenna number (1-14)
                     'mount' A string indicating the mount type for this antenna ('EQ' or 'AZEL')
                     'dra'   Measured RA offset (degrees) [npt]
                     'ddec'  Measured Dec offset (degrees) [npt]
                     'dxel'  XEL offset calculated from measured offsets (degrees) [npt]
                     'd_el'  EL offset calculated from measured offsets (degrees) [npt]
                     'ra'    Right ascension of antenna during measurement (radians) [npt]
                     'dec'   Declination of antenna during measurement (radians) [npt]
                     'az'    Azimuth of antenna during measurement (radians) [npt]
                     'el'    Elevation of antenna during measurement (radians) [npt]
                     'ha'    Hour angle of antenna during measurement (radians) [npt]

        Method:
          For each point, calculates refraction-corrected elevation el', where
              el' = el+(1/60.)*(0.0019279 + 1.02/tan((el + 10.3/(el+5.1))*dtor))   ## Skip refraction correction!
          then for an AZEL mount, does
              x = [1, cos(el'), sin(el'), cos(az)*sin(el'), sin(az)*sin(el'),  # Xel coordinates (angle on sky!)
                   0, 0, 0]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              x = [0, 0, 0, -sin(az), cos(az), 1, cos(el'), cot(el')]          # El coordinates
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
          while for an EQ mount, does
              x = [1, -cos(lat)*sin(ha)/cos(dec), tan(dec), -1./cos(dec),      # HA coordinates (not angle on sky!)
                   sin(ha)*tan(dec), -cos(ha)*tan(dec), 0, 0]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              x = [0, 0, 0, 0, cos(ha), sin(ha), 1,                            # Declination coordinates
                   -cos(lat)*cos(ha)*sin(dec) + sin(lat)*cos(dec)]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              
          Then
              p = np.linalg.solve(a,b)
              
          where p is the array of pointing parameters with different meanings according to mount type
              AZEL:                                          EQ:
              -------------------------------------------    -----------------------------------------
              P1 = Azimuth collimation error                 = Hour angle encoder offset
              P2 = Azimuth encoder offset                    = Hour angle sag
              P3 = Elevation axis skew angle                 = Axis skew
              P4 = -phi*sin(Ka), where phi = azimuth axis    = Collimation error
                   tilt angle and Ka = angle defining 
                   direction of tilt
              P5 = phi*cos(Ka)                               = Tilt out
              P6 = Not used                                  = Tilt over
              P7 = Elevation encoder offset and              = Declination encoder offset
                   collimation error                           and collimation error
              P8 = Gravitational deflection coefficient      = Declination sag
              P9 = Residual refraction coefficient           = Not used
    '''
    from numpy import cos, sin, tan  # Avoids all the np. etc.
    mount = indict['mount']
    npt = len(indict['az'])
    nparm = 8
    a = np.asmatrix(np.zeros((nparm,nparm))) # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)                      # Vector to create from measurements
    # Do accumulation over sources
    for i in range(npt):
        if mount == 'AZEL':
            el = indict['el'][i]
            az = indict['az'][i]
            dxel = indict['dxel'][i]
            d_el = indict['d_el'][i]
            elp = el   # Skip refraction correction
            # Add these pointing results to A and B
            x = np.array([1., np.cos(el), sin(el), cos(az)*sin(el), sin(az)*sin(el), 0., 0., 0.])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += dxel * x
            x = np.array([0., 0., 0., -sin(az), cos(az), 1., cos(el), 1./tan(el)])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += d_el * x
        elif mount == 'EQ':
            ha = indict['ha'][i]
            dec = indict['dec'][i]
            dha = -indict['dra'][i]   # Change sign for RA -> HA
            ddec = indict['ddec'][i]
            if indict['ant'] == 14:
                # Special case for broken pointing model in Ant 14 [tan(ha) is calculated instead of tan(dec)]
                x = np.array([1, -cos(lat)*sin(ha)/cos(dec), tan(ha), -1./cos(dec), 
                              sin(ha)*tan(ha), -cos(ha)*tan(ha), 0., 0.])
            else:
                x = np.array([1, -cos(lat)*sin(ha)/cos(dec), tan(dec), -1./cos(dec), 
                              sin(ha)*tan(dec), -cos(ha)*tan(dec), 0., 0.])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += dha * x
            x = np.array([0., 0., 0., 0., cos(ha), sin(ha), 1., -cos(lat)*cos(ha)*sin(dec) + sin(lat)*cos(dec)])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += ddec * x
            
    # Calculate the solution (p are best-fit parameters, P1-P9)
    p = np.linalg.solve(a,b)
    # Add parameters to indict and return    
    indict.update({'params':p}) #, 'el_r':np.array(el_r)}) # Skip refraction correction
    return indict
    
def checkfit(indict):
    ''' Check the fit from the dictionary returned by mntcal().
        
        Pointing model depends on mount type.
        AZEL mount type:  Predicted pointing offsets are

          XELO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
          ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)

        where XELO, ELO are small, measured pointing offsets in "azimuth" (i.e. cross-elevation) 
        and elevation angles on the sky, at azimuth, AZ, and elevation EL.

        EQ mount type: Predicted pointing offsets are
          HAO  = P1 - P2 cos(LAT)*sin(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) 
               + P5 sin(HA)*tan(DEC) - P6 cos(HA)*tan(DEC)
          DECO = P5 cos(HA) + P6 sin(HA) + P7 - P8 [cos(LAT)*cos(HA)*sin(DEC) - sin(LAT)*cos(DEC)]

        where HAO, DECO are small, measured pointing offsets in hour-angle and elevation
        (HA not angle on the sky), at hour-angle, HA, and declination EL.

    '''
    import matplotlib
    import matplotlib.pylab as plt
    import os
    from numpy import cos, sin, tan  # Avoids all the np. etc.
    mount = indict['mount']
    npt = len(indict['az'])
    nparm = 8
    p = indict['params']
    old_p = indict['params_old']
    fit_x = []   # Fit of horizontal offset (HA or AZ)
    fit_y = []   # Fit of vertical offset (DEC or EL)
    # Do fit to all of the measurements using the solution for P
    if mount == 'AZEL':
        # AZO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
        # ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)
        # Note: since P6 is not used, and indexes below are 0-based, these become
        # ELO = P[5] + P[6] cos(EL') + P[7] cot(EL') - P[3] sin(AZ) + P[4] cos(AZ)
        for i in range(npt):
            #elp  = indict['el_r'][i]
            elp  = indict['el'][i]
            az   = indict['az'][i]

            xel = p[0] + p[1]*cos(elp) + p[2]*sin(elp) + p[3]*cos(az)*sin(elp) + p[4]*sin(az)*sin(elp)
            fit_x.append(xel)
            el = p[5] + p[6]*cos(elp) + p[7]/tan(elp) - p[3]*sin(az) + p[4]*cos(az)
            fit_y.append(el)

        # Calculate the difference between the measured offsets and the fitted ones
        fit_x = np.array(fit_x)
        fit_y = np.array(fit_y)
        x_diff = indict['dxel'] - fit_x
        y_diff = indict['d_el'] - fit_y
    elif mount == 'EQ':
        # HAO  = P1 - P2 cos(LAT)*sin(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) 
        #      - P6 cos(HA)*tan(DEC)
        # DECO = P5 cos(HA) + P6 sin(HA) + P7 - P8 [cos(LAT)*cos(HA)*sin(DEC) - sin(LAT)*cos(DEC)]
        #      + P10 DEC
        for i in range(npt):
            dec  = indict['dec'][i]
            ha   = indict['ha'][i]
            if indict['ant'] == 14:
                # Special case for broken pointing model in Ant 14 [tan(ha) is calculated instead of tan(dec)]
                dh = p[0] - p[1]*cos(lat)*sin(ha)/cos(dec) + p[2]*tan(ha) - p[3]/cos(dec) \
                          + p[4]*sin(ha)*tan(ha) - p[5]*cos(ha)*tan(ha)
            else:
                dh = p[0] - p[1]*cos(lat)*sin(ha)/cos(dec) + p[2]*tan(dec) - p[3]/cos(dec) \
                          + p[4]*sin(ha)*tan(dec) - p[5]*cos(ha)*tan(dec)
            fit_x.append(dh)
            dd = p[4]*cos(ha) + p[5]*sin(ha) + p[6] - p[7]*(cos(lat)*cos(ha)*sin(dec) - sin(lat)*cos(dec))
            fit_y.append(dd)

        # Calculate the difference between the measured offsets and the fitted ones
        fit_x = np.array(fit_x)
        fit_y = np.array(fit_y)
        x_diff = -indict['dra'] - fit_x  # Change sign for RA -> HA
        y_diff = indict['ddec'] - fit_y

    # Print results
    print 'Solution for antenna',indict['ant'],'mount type:',indict['mount']
    print '                  Times \     AZ    EL     HA    DEC       MEASURED   FITTED  DIFFERENCE (deg)'
    print ' '
    if mount == 'AZEL':
        for i in range(npt):
            el  = indict['el'][i] / dtor
            az   = indict['az'][i] / dtor
            ha   = indict['ha'][i] / dtor
            dec  = indict['dec'][i] / dtor
            dxel = indict['dxel'][i]
            d_el = indict['d_el'][i]
            times = indict['times'][i]
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dXEL = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,dxel,fit_x[i],x_diff[i])
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dEL  = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,d_el,fit_y[i],y_diff[i])

        # Calculate an appropriate residual
        rmsum = np.concatenate((x_diff, y_diff)).std()
        origsum = np.sqrt((np.concatenate((indict['dxel'],indict['d_el']))**2).sum()/(npt*2 - nparm))

        # Add in zero for 6th parameter
        p = np.insert(p,5,0.0)
        # Print residual and solution
        print ' '
        print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
        print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'
    elif mount == 'EQ':   
        for i in range(npt):
            el  = indict['el'][i] / dtor
            az   = indict['az'][i] / dtor
            ha  = indict['ha'][i] / dtor
            dec   = indict['dec'][i] / dtor
            dha = -indict['dra'][i]
            ddec = indict['ddec'][i]
            times = indict['times'][i]
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dHA  = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,dha,fit_x[i],x_diff[i])
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dDec = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,ddec,fit_y[i],y_diff[i])

        # Calculate an appropriate residual
        rmsum = np.concatenate((x_diff, y_diff)).std()
        origsum = np.sqrt((np.concatenate((-indict['dra'],indict['ddec']))**2).sum()/(npt*2 - nparm))

        # Add in zero for 9th parameter
        p = np.insert(p,8,0.0)
        # Print residual and solution
        print ' '
        print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
        print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'

    # Add update to existing pointing parameters and provide new parameters
    print '{:7d}{:7d}{:7d}{:7d}{:7d}{:7d}{:7d}{:7d}{:7d}'.format(*(old_p + (p*10000).astype(int)))+' \ Updated pointing parameters'
    print ''
    indict.update({'params_new':old_p + (p*10000).astype(int)})
    '''
    #    paz = Azimuth list for plotting
    #    pel = Elevation list for plotting
    #    azpo = Az pointing offset for plotting
    #    elpo = El pointing offset for plotting
    #    azfit = Fit for Azimuth for plotting
    #    elfit = Fit for Elevation for plotting
        # Azimuth -- get index for sort by ascending order
        print 'start of plot'
    '''
    paz = indict['az'] / dtor
    pel = indict['el'] / dtor
    ind = paz.argsort()
    if mount == 'AZEL':
        azpo = indict['dxel']
        elpo = indict['d_el']
    elif mount == 'EQ':
        azpo = -indict['dra']
        elpo = indict['ddec']
    drange = np.sqrt((azpo.max() - azpo.min())**2 + (elpo.max() - elpo.min())**2)
    matplotlib.rcParams.update({'font.size':12})
    f = plt.figure()
    f.suptitle(indict['times'][0].iso[:16]+'-'+indict['times'][-1].iso[11:16]+'    Ant '+str(indict['ant']))
    plt.subplot(221)
    if drange > 0.25:
        plt.axis([30,330,-1.5,1.5])
    else:
        plt.axis([30,330,-0.25,0.25])
    plt.xlabel('Azimuth [deg]')
    if mount == 'AZEL':
        plt.ylabel('Azimuth Offset [deg]')
        plt.title('Azimuth Offsets and Fit')
    elif mount == 'EQ':
        plt.ylabel('HA Offset [deg]')
        plt.title('HA Offsets and Fit')
    plt.plot(paz[ind],azpo[ind],'o')
    plt.plot(paz[ind],fit_x[ind])

    plt.subplot(222)
    if drange > 0.25:
        plt.axis([30,330,-1.5,1.5])
    else:
        plt.axis([30,330,-0.25,0.25])
    plt.xlabel('Azimuth [deg]')
    if mount == 'AZEL':
        plt.ylabel('Elevation Offset [deg]')
        plt.title('Elevation Offsets and Fit')
    elif mount == 'EQ':
        plt.ylabel('Declination Offset [deg]')
        plt.title('Declination Offsets and Fit')
    plt.plot(paz[ind],elpo[ind],'o')
    plt.plot(paz[ind],fit_y[ind])

    plt.subplot(223)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([30,300,10,88])
    plt.xlabel('Azimuth [deg]')
    plt.ylabel('Elevation [deg]')
    plt.title('Sky Coverage')
    plt.plot(paz[ind],pel[ind],'+')

    plt.subplot(224)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-0.5,0.5,-0.5,0.5])
    if mount == 'AZEL':
        plt.xlabel('Azimuth Offset [deg]')
        plt.ylabel('Elevation Offset [deg]')
    elif mount == 'EQ':
        plt.xlabel('HA Offset [deg]')
        plt.ylabel('Declination Offset [deg]')
    plt.text(0.0,0.4,'Pointing vs. Solar Disk',ha='center')
    plt.text(-0.45,-0.38,'$RMS_{orig}$',color='C1',ha='left')
    plt.text(0.45,-0.38,'$RMS_{corr}$',color='C3',ha='right')
    plt.text(-0.45,-0.45,'{:5.3f}'.format(origsum)+'$^o$',color='C1',ha='left')
    plt.text(0.45,-0.45,'{:5.3f}'.format(rmsum)+'$^o$',color='C3',ha='right')
    th = np.linspace(0,2*np.pi,100)
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),azpo,elpo,'+')
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),x_diff,y_diff,'.')

    # Set plot size
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    # Save figure
    filename = indict['filename']
    path = os.path.dirname(filename)
    if path == '': path = '.'
    fname = os.path.basename(filename)[:8]+'-ant'+str(indict['ant'])+'-pnt.png'
    plt.savefig(path+os.sep+fname)
    plt.close()
    return indict

def multi_mountcal(filename, ant_str=None):
    ''' Process an entire set of calpnt data
    '''
    import coord_conv as cc
    from util import ant_str2list
    
    outdict = rd_calpnt(filename)
    indict_list = []
    if ant_str:
        antlist = ant_str2list(ant_str)+1
    else:
        antlist = outdict['antlist']
    for ant in antlist:
        #i = ant - 1
        i, = np.where(np.array(outdict['antlist']) == ant)[0]  # Index in outdict for specified ant
        if ant in [9, 10, 11, 13, 14]:
            # This is an equatorial mount
            mount = 'EQ'
        elif ant in [1, 2, 3, 4, 5, 6, 7, 8, 12]:
            # This is an azimuth-elevation mount
            mount = 'AZEL'
        # Reprocess coordinates, plus remove any -99 points
        good, = np.where(np.logical_and(outdict['dra'][i] > -90, outdict['ddec'][i] > -90))
        npt = len(good)
        az = np.zeros(npt)
        el = np.zeros(npt)
        dxel = np.zeros(npt)
        d_el = np.zeros(npt)
        dra = outdict['dra'][i,good]
        ddec = outdict['ddec'][i,good]
        ra = outdict['ra'][good]
        dec = outdict['dec'][good]
        ha = outdict['ha'][good]
        times = outdict['time'][good]
        params_old = outdict['params_old'][:,i]
        # el_r = []
        for k in range(npt):
            # Convert RA, Dec to Az, El, and dRA, dDec to dxel and d_el
            azk, elk = cc.radec2azel(ra[k], dec[k], times[k])
            # Apply refraction correction (assumes elevation in degrees)
            # This is commented out, pending determination of whether refraction is already
            # accounted for in the measurements (I think it is...)
            #elk /= dtor  # Convert to degrees
            #elp = elk+(1/60.)*(0.0019279 + 1.02/tan((elk + 10.3/(el+5.1))*dtor))
            #elp *= dtor # Convert back to radians
            #el_r.append(elp)
            # Adjust coordinates (only needed for AZEL, but do for EQ, too for consistency)
            dxelk, d_elk = cc.dradec2dazel(ra[k],dec[k],times[k],dra[k]*dtor,ddec[k]*dtor)
            az[k] = azk
            el[k] = elk
            dxel[k] = dxelk / dtor
            d_el[k] = d_elk / dtor
        indict = {'filename':outdict['filename'],'ant':ant, 'dra':dra, 'ddec':ddec, 'dxel': dxel, 'd_el':d_el, 'ra':ra, 
                  'dec':dec, 'ha':ha, 'az':az, 'el':el, 'mount':mount, 'times':times, 
                  'params_old':params_old}
        # indict.update({'el_r':el_r})    # Refraction-corrected elevation
        indict = mntcal(indict)
        indict = checkfit(indict)
        indict_list.append(indict.copy())
    return indict_list   # Temporary
    
def params2ants(indict_list):
    ''' Send params_new values from indict_list (output of multi_mountcal) to antennas
    
        indict_list[i]['ant']           gives the antenna to send information to
        indict_list[i]['params_new']    contains the parameters to send
        
    '''
    import calibration as cal
    import stateframe as stf
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    for indict in indict_list:
        ant = str(indict['ant'])
        for i in range(9):
            cmd = 'pointingcoefficient'+str(i+1)+' '+str(indict['params_new'][i])+' ant'+ant
            cal.send_cmds([cmd],acc)
        print 'Coefficients for ant',ant,'sent.'
        if int(ant) in [1,2,3,4,5,6,7,8,12,14]:
            print 'Antenna',ant,'controller must be rebooted before parameters take effect.'
