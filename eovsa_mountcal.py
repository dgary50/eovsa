from util import Time
import numpy as np
from eovsa_lst import eovsa_lst

def rd_calpnt(filename):
    ''' Read and return contents of output of CALPNT or CALPNT2M observation
    
    '''
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
    return {'filename':filename, 'source':source, 'time':times, 'ra':ra, 
            'dec':dec, 'ha':ha, 'antlist':antlist, 'dra':dra, 'ddec':ddec}
            
def mntcal(indict):
    ''' Given an input dictionary of pointing coordinates for a single AZEL antenna
        calculate pointing parameters.
        
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
              x = [1, cos(el'), sin(el'), cos(az)*sin(el'), sin(az)*sin(el'), 0, 0, 0]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              x = [0, 0, 0, -sin(az), cos(az), 1, cos(el'), cot(el')]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
          while for an EQ mount, does
              x = [1, -cos(lat)*sin(ha)/cos(dec), tan(dec), -1./cos(dec), sin(ha)*tan(dec), -cos(ha)*tan(dec), 0, 0]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              x = [0, 0, 0, 0, cos(ha), sin(ha), 1, -cos(lat)*cos(ha)*sin(dec) + sin(lat)*cos(dec)]
              a += asmatrix(x).T * asmatrix(x)
              b += dxel * x
              
          Then
              p = np.linalg.solve(a,b)
              
          where p is the array of pointing parameters
              P1 = Azimuth collimation error
              P2 = Azimuth encoder offset
              P3 = Elevation axis skew angle
              P4 = -phi*sin(Ka), where phi = azimuth axis tilt angle and Ka = angle defining direction of tilt
              P5 =  phi*cos(Ka)
              P6 (not used)
              P7 = Elevation encoder offset + collimation error
              P8 = Gravitational deflection coefficient
              P9 = Residual refraction coefficient        
    '''
    from numpy import cos, sin, tan  # Avoids all the np. etc.
    dtor = np.pi/180.
    npt = len(indict['az'])
    nparm = 8
    a = np.asmatrix(np.zeros((nparm,nparm))) # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)                      # Vector to create from measurements
    el_r = []    # List of refraction-corrected elevations
    # Do accumulations
    for i in range(npt):
        
        # Apply refraction correction (assumes elevation in degrees)
        #el /= dtor  # Convert to degrees
        #elp = el+(1/60.)*(0.0019279 + 1.02/tan((el + 10.3/(el+5.1))*dtor))
        #elp *= dtor # Convert to radians
        #el_r.append(elp)
        if mount == 'AZEL':
            el = indict['el'][i]
            az = indict['az'][i]
            dxel = indict['dxel'][i]
            d_el = indict['d_el'][i]
            elp = el   # Skip refraction correction
            # Add these pointing results to A and B
            x = np.array([1., np.cos(elp), sin(elp), cos(az)*sin(elp), sin(az)*sin(elp), 0., 0., 0.])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += dxel * x
            x = np.array([0., 0., 0., -sin(az), cos(az), 1., cos(elp), 1./tan(elp)])
            a += np.asmatrix(x).T * np.asmatrix(x)
            b += d_el * x
        elif mount == 'EQ':
            ha = indict['ha'][i]
            dec = indict['dec'][i]
            dha = -indict['dra'][i]   # Change sign for RA -> HA
            ddec = indict['ddec'][i]
            x = np.array([1, -cos(lat)*sin(ha)/cos(dec), tan(dec), -1./cos(dec), sin(ha)*tan(dec), -cos(ha)*tan(dec), 0., 0.])
            a += asmatrix(x).T * asmatrix(x)
            b += dha * x
            x = np.array([0., 0., 0., 0., cos(ha), sin(ha), 1., -cos(lat)*cos(ha)*sin(dec) + sin(lat)*cos(dec)])
            a += asmatrix(x).T * asmatrix(x)
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

          AZO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
          ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)

        where AZO, ELO are small, measured pointing offsets in "azimuth" and elevation
        angles on the sky, at azimuth, AZ, and elevation EL.

        EQ mount type: Predicted pointing offsets are
          HAO  = P1 - P2 cos(LAT)*sin(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) - P6 cos(HA)*tan(DEC)
          DECO = P5 cos(HA) + P6 sin(HA) + P7 - P8 [cos(LAT)*cos(HA)*sin(DEC) - sin(LAT)*cos(DEC)]

        where HAO, DECO are small, measured pointing offsets in hour-angle and elevation
        angles on the sky, at hour-angle, HA, and declination EL.

    '''
    import matplotlib
    import matplotlib.pylab as plt
    from numpy import cos, sin, tan  # Avoids all the np. etc.
    dtor = np.pi/180.
    npt = len(indict['az'])
    nparm = 8
    p = indict['params']
    fit_x = []   # Fit of horizontal offset (HA or AZ)
    fit_y = []   # Fit of vertical offset (DEC or EL)
    # Do fit to all of the measurements using the solution for P
    if mount == 'AZEL':
        # AZO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
        # ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)
        for i in range(npt):
            elp  = indict['el_r'][i]
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
            
            dh = p[0] - p[1]*cos(lat)*sin(ha[i])/cos(dec[i]) + p[2]*tan(dec[i]) - p[3]/cos(dec[i]) \
                      + p[4]*sin(ha[i])*tan(dec[i]) - p[5]*cos(ha[i])*tan(dec[i])
            fit_x.append(dh)
            dd = p[4]*cos(ha[i]) + p[5]*sin(ha[i]) + p[6] - p[7]*(cos(lat)*cos(ha[i])*sin(dec[i]) - sin(lat)*cos(dec[i]))
            fit_y.append(dd)

        # Calculate the difference between the measured offsets and the fitted ones
        fit_x = np.array(fit_x)
        fit_y = np.array(fit_y)
        x_diff = -indict['dra'] - fit_x  # Change sign for RA -> HA
        y_diff = indict['ddec'] - fit_y

    # Print results
    print ' Solution for antenna',indict['ant'],'mount type:',indict['mount']
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
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dXEL = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,dxel,fit_x[i],xel_diff[i])
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dEL  = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,d_el,fit_y[i],el_diff[i])

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
            dra = indict['dra'][i]
            ddec = indict['ddec'][i]
            times = indict['times'][i]
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dHA = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,dha,fit_x[i],x_diff[i])
            print times, "\ {:6.1f}{:6.1f} {:6.1f}{:6.1f}   dDec  = {:7.3f}  {:7.3f}  {:7.3f}".format(az,el,ha,dec,ddec,fit_y[i],y_diff[i])

        # Calculate an appropriate residual
        rmsum = np.concatenate((x_diff, y_diff)).std()
        origsum = np.sqrt((np.concatenate((indict['dha'],indict['ddec']))**2).sum()/(npt*2 - nparm))

        # Add in zero for 6th parameter
        p = np.insert(p,5,0.0)
        # Print residual and solution
        print ' '
        print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
        print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'

    '''    
    print ''.join('{:7d}'.format(int(int(aligntab[k])+v*10000)) for k,v in enumerate(p)),' ALIGNPARM \ Updated pointing parameters'

#    paz = Azimuth list for plotting
#    pel = Elevation list for plotting
#    azpo = Az pointing offset for plotting
#    elpo = El pointing offset for plotting
#    azfit = Fit for Azimuth for plotting
#    elfit = Fit for Elevation for plotting
    # Azimuth -- get index for sort by ascending order
    print 'start of plot'
    '''
    # paz = indict['az'] / dtor
    # pel = indict['el'] / dtor
    # ind = paz.argsort()
    # azpo = indict['dxel']
    # elpo = indict['d_el']
    # print 'setting drange'
    # drange = np.sqrt((azpo.max() - azpo.min())**2 + (elpo.max() - elpo.min())**2)
    # matplotlib.rcParams.update({'font.size':12})
    # plt.figure()
    # plt.subplot(221)
    # if drange > 0.25:
        # plt.axis([30,330,-1.5,1.5])
    # else:
        # plt.axis([30,330,-0.25,0.25])
    # plt.xlabel('Azimuth [deg]')
    # plt.ylabel('Azimuth Offset [deg]')
    # plt.title('Azimuth Offsets and Fit')
    # plt.plot(paz[ind],azpo[ind],'o')
    # plt.plot(paz[ind],fit_x[ind])

    # plt.subplot(222)
    # if drange > 0.25:
        # plt.axis([30,330,-1.5,1.5])
    # else:
        # plt.axis([30,330,-0.25,0.25])
    # plt.xlabel('Azimuth [deg]')
    # plt.ylabel('Elevation Offset [deg]')
    # plt.title('Elevation Offsets and Fit')
    # plt.plot(paz[ind],elpo[ind],'o')
    # plt.plot(paz[ind],fit_y[ind])

    # plt.subplot(223)
    # plt.axis('equal')
    # plt.axis('scaled')
    # plt.axis([30,300,10,88])
    # plt.xlabel('Azimuth [deg]')
    # plt.ylabel('Elevation [deg]')
    # plt.title('Sky Coverage')
    # plt.plot(paz[ind],pel[ind],'+')

    # plt.subplot(224)
    # plt.axis('equal')
    # plt.axis('scaled')
    # plt.axis([-0.5,0.5,-0.5,0.5])
    # plt.xlabel('Azimuth Offset [deg]')
    # plt.ylabel('Elevation Offset [deg]')
    # plt.title('Pointing Relative to Solar Disk')
    # th = np.linspace(0,2*np.pi,100)
    # plt.plot(0.25*np.cos(th),0.25*np.sin(th),azpo,elpo,'+')
    # # Set plot size
    # #fig = plt.gcf()
    # #fig.set_size_inches(10.0,8.0)
    # #tok = filename.split('/')
    # #stem = tok[len(tok)-1].split('.')[0]
    # #plt.savefig(stem+'.pdf')
    # #plt.close()

def multi_mountcal(filename):
    ''' Process an entire set of calpnt data
    '''
    import coord_conv as cc
    outdict = rd_calpnt(filename)
    indict_list = []
    for i, ant in enumerate(outdict['antlist']):
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
        for k in range(npt):
            # Convert RA, Dec to Az, El, and dRA, dDec to dxel and d_el
            azk, elk = cc.radec2azel(ra[k], dec[k], times[k])
            # Adjust coordinates 
            dxelk, d_elk = cc.dradec2dazel(ra[k],dec[k],times[k],dra[k],ddec[k])
            az[k] = azk
            el[k] = elk
            dxel[k] = dxelk
            d_el[k] = d_elk
        indict = {'ant':ant, 'dra':dra, 'ddec':ddec, 'dxel': dxel, 'd_el':d_el, 'ra':ra, 'dec':dec, 'ha':ha, 'az':az, 'el':el, 'mount':mount, 'times':times}
        indict = mntcal(indict)
        indict_list.append(indict.copy())
        print 'Results for Antenna',ant
        checkfit(indict)
    return indict_list   # Temporary
                
