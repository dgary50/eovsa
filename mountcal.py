# History
#  2017-04-14  DG
#    Flipped the sign of p8 term in anta_mountcal() to reflect actual usage
#    in the Ant A controller.  I also added a temporary anta_broken_mountcal()
#    routine to fits to the "broken" pointing model actually in use in Ant A.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from readbsc import starobs2dxeldel
#+
# NAME:
#     KSRBL_MOUNTCAL
# PURPOSE:
#     Version of MOUNTCAL for KSRBL (Patriot) dish.  Based on OVSA Mountcal.
# CATEGORY:
#     KSRBL ANTENNA CALIBRATION
# CALLING SEQUENCE:
#     ksrbl_mountcal,filename[,/star]
# INPUTS:
#     filename  Name of file to read, e.g. starsolutions-2013-07-19.txt, containing 
#                  the measurement of stars, which are converted to offsets 
#                  internally by starobs2dxeldel().
# OPTIONAL (KEYWORD) INPUT PARAMETERS:
# ROUTINES CALLED:
# OUTPUTS:
# COMMENTS:
#     Algorithm is based on the following geometrical formulae
#        (from Patriot Documentation):
#
#     First, convert elevation of measurement to EL', corrected for refraction:
#
#                    1                      1.02
#        EL' = EL + --- ( 0.0019279 + ------------------)
#                   60.                          10.3
#                                     tan(EL + --------)
#                                              EL + 5.1
#
#     Then, invert system of equations:
#
#        AZO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
#        ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)
#
#     where AZO, ELO are small, measured pointing offsets in "azimuth" and elevation
#     angles on the sky, at azimuth, AZ, and elevation EL.
#
#     The demanded azimuth and elevation will be:
#
#                           AZO
#        AZ_demand = AZ + --------
#                         cos(EL')
#
#        EL_demand = EL' + ELO
#
#     P1 = Azimuth collimation error
#     P2 = Azimuth encoder offset
#     P3 = Elevation axis skew angle
#     P4 = -phi*sin(Ka), where phi = azimuth axis tilt angle and Ka = angle defining direction of tilt
#     P5 =  phi*cos(Ka)
#     P6 (not used)
#     P7 = Elevation encoder offset + collimation error
#     P8 = Gravitational deflection coefficient
#     P9 = Residual refraction coefficient
#
#     The input data is a set of measured values of AZO or ELO at different (AZ, EL).
#     The task is to find the best fit values of P1, P2, P3, P4, P5, P7, P8, P9.
#     This is equivalent to solving a set of 8 linear equations in 8 unknowns.
#
#     Let the i'th of NPT data points be:
#
#       Y[i]    =               AZO             or            ELO
#
#     Define:
#
#       X[1,i]  =               1               or              0
#       X[2,i]  =               cos(EL')        or              0
#       X[3,i]  =               sin(EL')        or              0
#       X[4,i]  =               cos(AZ)sin(EL') or       -sin(AZ)
#       X[5,i]  =               sin(AZ)sin(EL') or        cos(AZ)
#       X[6,i]  =               0               or              1
#       X[7,i]  =               0               or       cos(EL')
#       X[8,i]  =               0               or       cot(EL')
#
#     Then define the sums, A[j,k] = SUM over i  ( X[j,i]*X[k,i] )
#                           B[k]   = SUM over i  ( Y[i]  *X[k,i] )
#
#     Then the set of equations to be solved are:
#
#               SUM over j ( A[k,j] * Pj ) = B[k]               k = 1,2,3,4,5,6,7,8
#
# SIDE EFFECTS:
#     Writes results to the screen--not too useful
# RESTRICTIONS:
# MODIFICATION HISTORY:
#     Written 21-Oct-2007 by Dale E. Gary (modified from OVSA Mountcal)
#     30-Jun-2008  DG
#       Added /STAR keyword switch, which inverts the sign of the calculated pointing
#       parameters to account for the fact that optical pointing gives the position of
#       the pointing center, rather than the usual position of a source relative to the
#       pointing center.
#     12-Jan-2013  DG
#       Converted to Python
#     19-Jul-2013  DG
#       Added internal call to starobs2dxeldel()
#-
def mountcal(filename=None,param_string=None,star=True):

    nparm = 8

    # Open the file containing the raw pointing data
    if filename is None:
        print 'Error: Must specify an input file name.'
        return
    if param_string is None:
        param_string = '0 0 0 0 0 0 0 0 0'
    starobs2dxeldel(filename)
    idx = filename.find('solutions')
    infile = filename[0:idx]+'reduction'+filename[idx+9:]

    f = open(infile,'r')

    # Read first line to get old alignment parameters
    # line = f.readline()
    aligntab = param_string.strip().split() # Read old alignment parameters into ALIGNTAB
    # Read rest of file into lines
    lines = f.readlines()
    f.close()
    nlines = len(lines)

##    !p.multi = [0,2,1,0,0]
##    !p.charsize = 1

##    saveresult = ' '

     # Set some initializations
    n = 0          # Number of "good" measurements
    npt = 0        # Total number of entries (some could be "flagged" bad)

    x = np.zeros((nparm,nparm))
    a = np.asmatrix(x)   # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)  # Array to create from measurements
    aaz = np.zeros(nlines*2)
    paz = []
    pel = []
    azpo = []
    elpo = []
    ael = np.zeros(nlines*2)
    apo = np.zeros(nlines*2)
    ahodo = []
    hodo = []
    name = []
    dtor = np.pi/180.

    for line in lines:
        nam = line[:10]
        time,ra0,dec0,ra1,dec1,dra,ddec,az,el,daz,delv = line[10:].strip().split()

        # Enter line contents into arrays for later use
        name.append(nam)
        ahodo.append('AZO')
        aaz[npt] = float(az)
        paz.append(float(az))
        el = float(el)
        pel.append(float(el))
        ael[npt] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
        apo[npt] = float(daz)
        azpo.append(float(daz))
        name.append(nam)
        ahodo.append('ELO')
        aaz[npt+1] = float(az)
        ael[npt+1] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
        apo[npt+1] = float(delv)
        elpo.append(float(delv))
        
        # Coordinate parameters are contained in variable X, which are different
        # for azimuth and elevation measurements
        if ahodo[npt] is 'AZO':
             x = np.array([1.,
                 np.cos(ael[npt]*dtor),
                 np.sin(ael[npt]*dtor),
                 np.cos(aaz[npt]*dtor)*np.sin(ael[npt]*dtor),
                 np.sin(aaz[npt]*dtor)*np.sin(ael[npt]*dtor),
                 0.,
                 0.,
                 0.])
             n += 1
             b += np.array(apo[npt]*x)
             xx = []
             for i in range(len(x)):
                 xx.append(x*x[i])
             a += np.asmatrix(xx)
        if ahodo[npt+1] is 'ELO':
             x = np.array([0.,
                  0.,
                  0.,
                 -np.sin(aaz[npt+1]*dtor),
                  np.cos(aaz[npt+1]*dtor),
                  1.,
                  np.cos(ael[npt+1]*dtor),
                  1./np.tan(ael[npt+1]*dtor)])
             n += 1
             b += np.array(apo[npt+1]*x)
             xx = []
             for i in range(len(x)):
                 xx.append(x*x[i])
             a += np.asmatrix(xx)

        npt += 2

    if npt == 0:
        print 'KSRBL_MOUNTCAL: File read error.'
        return

    # Solve equation for pointing parameters
    p = np.linalg.solve(a,b)

    fit = []
    azfit = []
    elfit = []
    # Do fit to all of the measurements using the solution for P
    # AZO = P1 + P2 cos(EL') + P3 sin(EL') + P4 cos(AZ)sin(EL') + P5 sin(AZ)sin(EL')
    # ELO = P7 + P8 cos(EL') + P9 cot(EL') - P4 sin(AZ)         + P5 cos(AZ)
    for i in range(npt):
        if ahodo[i] is 'AZO':
            az = p[0] + p[1]*np.cos(ael[i]*dtor) \
                      + p[2]*np.sin(ael[i]*dtor) \
                      + p[3]*np.cos(aaz[i]*dtor)*np.sin(ael[i]*dtor) \
                      + p[4]*np.sin(aaz[i]*dtor)*np.sin(ael[i]*dtor)
            fit.append(az)
            azfit.append(az)
        elif ahodo[i] is 'ELO':
            el = p[5] + p[6]*np.cos(ael[i]*dtor) \
                      + p[7]/np.tan(ael[i]*dtor) \
                      - p[3]*np.sin(aaz[i]*dtor) \
                      + p[4]*np.cos(aaz[i]*dtor)
            fit.append(el)
            elfit.append(el)

    # Calculate the difference between the measured offsets and the fitted ones
    diff = np.array(apo) - np.array(fit)

    # Print results
    print ' '
    print '\     AZ    EL       MEASURED   FITTED  DIFFERENCE (deg)'
    print ' '
    for i in range(npt):
        print "\ {:6.1f}{:6.1f}   {:s}= {:7.3f}  {:7.3f}  {:7.3f}".format(aaz[i],ael[i],ahodo[i],apo[i],fit[i],diff[i])

    # Calculate an appropriate residual
    rmsum = diff.std()
    origsum = (apo**2).sum()
    origsum = np.sqrt(origsum/(n-nparm))

    # If these are optical stellar measurements, change the sign of the corrections,
    # since they are positions of the center of the field, not the offset of the source
    # from the center.
    if star: p = -p

    # Add in zero for 6th parameter
    p = np.insert(p,5,0.0)
    # Print residual and solution
    print ' '
    print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
    print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'
    
    print ''.join('{:7d}'.format(int(int(aligntab[k])+v*10000)) for k,v in enumerate(p)),' ALIGNPARM \ Updated pointing parameters'

#    paz = Azimuth list for plotting
#    pel = Elevation list for plotting
#    azpo = Az pointing offset for plotting
#    elpo = El pointing offset for plotting
#    azfit = Fit for Azimuth for plotting
#    elfit = Fit for Elevation for plotting
    # Azimuth -- get index for sort by ascending order
    print 'start of plot'
    paz = np.array(paz)
    pel = np.array(pel)
    ind = paz.argsort()
    ind_el = pel.argsort()
    azpo = np.array(azpo)
    elpo = np.array(elpo)
    print 'setting drange'
    drange = np.sqrt((azpo.max() - azpo.min())**2 + (elpo.max() - elpo.min())**2)
    matplotlib.rcParams.update({'font.size':12})
    plt.subplot(221)
    if drange > 0.25:
        plt.axis([30,330,-1.5,1.5])
    else:
        plt.axis([30,330,-0.25,0.25])
    #plt.xlabel('Azimuth [deg]')
    plt.xlabel('Elevation [deg]')
    plt.ylabel('Azimuth Offset [deg]')
    plt.title('Azimuth Offsets and Fit')
    #plt.plot(paz[ind],np.array(azpo)[ind],'o')
    #plt.plot(paz[ind],np.array(azfit)[ind])
    plt.plot(pel[ind_el],np.array(azpo)[ind_el],'o')
    plt.plot(pel[ind_el],np.array(azfit)[ind_el])

    plt.subplot(222)
    if drange > 0.25:
        plt.axis([30,330,-1.5,1.5])
    else:
        plt.axis([30,330,-0.25,0.25])
    plt.xlabel('Azimuth [deg]')
    plt.ylabel('Elevation Offset [deg]')
    plt.title('Elevation Offsets and Fit')
    plt.plot(paz[ind],np.array(elpo)[ind],'o')
    plt.plot(paz[ind],np.array(elfit)[ind])

    plt.subplot(223)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([30,300,10,88])
    plt.xlabel('Azimuth [deg]')
    plt.ylabel('Elevation [deg]')
    plt.title('Sky Coverage')
    plt.plot(paz[ind],np.array(pel)[ind],'+')

    plt.subplot(224)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-0.25,0.25,-0.25,0.25])
    plt.xlabel('Azimuth Offset [deg]')
    plt.ylabel('Elevation Offset [deg]')
    plt.title('Pointing Relative to Solar Disk')
    th = np.linspace(0,2*np.pi,100)
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),np.array(azpo),np.array(elpo),'+',color='C0')
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),np.array(azpo)-np.array(azfit),np.array(elpo)-np.array(elfit),'+',color='C1')
    # Set plot size
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    tok = filename.split('/')
    stem = tok[len(tok)-1].split('.')[0]
    plt.savefig(stem+'.pdf')
    #plt.close()
    

#+
# NAME:
#     EQ_MOUNTCAL
# PURPOSE:
#     Version of MOUNTCAL for equatorial dishes.  Based on the above Mountcal.
# CATEGORY:
#     ANTENNA CALIBRATION
# CALLING SEQUENCE:
#     eq_mountcal,filename[,star=True][,stepsize=stepsize]
# INPUTS:
#     filename  Name of file to read, e.g. starsolutions-2015-10-20.txt, containing 
#                  the measurement of stars, which are converted to offsets 
#                  starobs2dhaddec().
# OPTIONAL (KEYWORD) INPUT PARAMETERS:
#     star   Boolean (default True) determines whether pointing information is
#              obtained from star pointing information, or radio pointing.
# ROUTINES CALLED:
# OUTPUTS:
# COMMENTS:
#     Algorithm is based on the following geometrical formulae
#        (from InterTronic Solutions Documentation):
#
#     First, convert elevation of measurement to EL', corrected for refraction:
#
#                    1                      1.02
#        EL' = EL + --- ( 0.0019279 + ------------------)
#                   60.                          10.3
#                                     tan(EL + --------)
#                                              EL + 5.1
#
#     Then, invert system of equations:
#
#        HAO = P1 - P2 cos(LAT)*cos(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) 
#                 - P6 cos(HA)*tan(DEC) + P9 HA
#        DECO = P5 cos(HA) + P6 sin(HA) + P7 + P8 [cos(LAT)*cos(HA)*sin(DEC) + sin(LAT)*cos(DEC)]
#                 + P10 DEC
#
#     where HAO, DECO are small, measured pointing offsets in hour angle (inverse sign 
#     of differences in RA) and declination angles on the sky, at hour angle HA and declination 
#     DEC, for location latitude LAT.
#
#     The demanded azimuth and elevation will be:
#
#                           AZO
#        AZ_demand = AZ + --------
#                         cos(EL')
#
#        EL_demand = EL' + ELO
#
#     P1 = Hour angle offset
#     P2 = Hour angle sag
#     P3 = Axis skew angle
#     P4 = Collimation error
#     P5 =  Tilt out
#     P6 = Tilt over
#     P7 = Declination encoder offset + collimation error
#     P8 = Declination sag
#     P9 = HA step size factor (> 0 => stepsize is too small)   \  Note that P9 and P10 should be equal
#     P10 = Dec step size factor (> 0 => stepsize is too small) /
#
#     The input data are a set of measured values of RAO (-HAO) or DECO at different (HA, DEC).
#     The task is to find the best fit values of P1, P2, P3, P4, P5, P6, P7, P8, P9, P10.
#     This is equivalent to solving a set of 10 linear equations in 10 unknowns.  If stepsize
#     is not given, then the solution is for P1-P8 (8 equations in 8 unknowns), which is the
#     preferred solution assuming that the stepsize is known.
#
#     Let the i'th of NPT data points be:
#
#       Y[i]    =               HAO             or            DECO
#
#     Define:
#
#       X[1,i]  =               1               or              0
#       X[2,i]  =    -cos(LAT)sin(HA)sec(DEC)   or              0
#       X[3,i]  =            tan(DEC)           or              0
#       X[4,i]  =           -sec(DEC)           or              0
#       X[5,i]  =        sin(HA)tan(DEC)        or            cos(HA)
#       X[6,i]  =       -cos(HA)tan(DEC)        or            sin(HA)
#       X[7,i]  =               0               or              1
#       X[8,i]  =               0               or  -[cos(LAT)cos(HA)sin(DEC) - sin(LAT)cos(DEC)]
#       X[9,i]  =              HA               or              0
#       X[10,i] =               0               or             DEC-LAT
#
#  NB: The DECO X[8,i] term has signs that differ from the 27-m Antenna Controller document
#      provided by InterTronic Solutions.  After an exchange with Mark Godwin, it appears that
#      the documentation is wrong, and the code implemented in the controller is right, with
#      signs as above.
#
#     Then define the sums, A[j,k] = SUM over i  ( X[j,i]*X[k,i] )
#                           B[k]   = SUM over i  ( Y[i]  *X[k,i] )
#
#     Then the set of equations to be solved are:
#
#               SUM over j ( A[k,j] * Pj ) = B[k]               k = 1,2,3,4,5,6,7,8,9,10
#
# MODIFICATION HISTORY:
#     Written 20-Oct-2015 by Dale E. Gary (modified from EOVSA Mountcal)
#-

def eq_mountcal(filename=None, param_string=None, star=True, stepsize=None, sfac=0.0):

    lat = 37.233170       # OVSA Latitude (degrees)
    nparm = 8
    # Open the file containing the raw pointing data
    if filename is None:
        print 'Error: Must specify an input file name.'
        return
    if stepsize != None:
        nparm = 10
    if param_string is None:
        if nparm == 10:
            param_string = '0 0 0 0 0 0 0 0 0 0'
        else:
            param_string = '0 0 0 0 0 0 0 0'
    starobs2dxeldel(filename,hadec=True)
    idx = filename.find('solutions')
    infile = filename[0:idx]+'reduction'+filename[idx+9:]

    f = open(infile,'r')

    # Read first line to get old alignment parameters
    # line = f.readline()
    aligntab = param_string.strip().split() # Read old alignment parameters into ALIGNTAB
    # Read rest of file into lines
    lines = f.readlines()
    f.close()
    nlines = len(lines)

##    !p.multi = [0,2,1,0,0]
##    !p.charsize = 1

##    saveresult = ' '

     # Set some initializations
    n = 0          # Number of "good" measurements
    npt = 0        # Total number of entries (some could be "flagged" bad)

    x = np.zeros((nparm,nparm))
    a = np.asmatrix(x)   # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)  # Array to create from measurements
    aha = np.zeros(nlines*2)
    pha = []
    pdec = []
    hapo = []
    decpo = []
    adec = np.zeros(nlines*2)
    apo = np.zeros(nlines*2)
    ahodo = []
    hodo = []
    name = []
    dtor = np.pi/180.

    for line in lines:
        nam = line[:10]
        time,ra0,dec0,ra1,dec1,dra,ddec,ha,dec,daz,delv = line[10:].strip().split()

        # If the stepsize factor is given as an argument, it is applied here.  What
        # this means is that if the stepsize were adjusted by sfac, the errors dra
        # and ddec would have been these adjusted values.
        dra = float(dra) - sfac*float(ha)
        ddec = float(ddec) + sfac*(float(dec)-lat)

        # Enter line contents into arrays for later use
        name.append(nam)
        ahodo.append('HAO')
        aha[npt] = float(ha)
        adec[npt] = float(dec)
        pha.append(float(ha))
        pdec.append(float(dec))
        # ael[npt] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
        apo[npt] = -dra     # Change sign for RA -> HA
        hapo.append(-dra)

        name.append(nam)
        ahodo.append('DECO')
        aha[npt+1] = float(ha)
        adec[npt+1] = float(dec)
        #ael[npt+1] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
        apo[npt+1] = ddec
        decpo.append(ddec)
        
        # Coordinate parameters are contained in variable X, which are different
        # for azimuth and elevation measurements
        if ahodo[npt] is 'HAO':
             x = [1.,
                 -np.cos(lat*dtor)*np.sin(aha[npt]*dtor)/np.cos(adec[npt]*dtor),
                 np.tan(adec[npt]*dtor),
                 -1./np.cos(adec[npt]*dtor),
                 np.sin(aha[npt]*dtor)*np.tan(adec[npt]*dtor),
                 -np.cos(aha[npt]*dtor)*np.tan(adec[npt]*dtor),
                 0.,
                 0.]
             if nparm == 10:
                 x.append(aha[npt])
                 x.append(0.)
             x = np.array(x)
             n += 1
             b += np.array(apo[npt]*x)
             xx = []
             for i in range(len(x)):
                 xx.append(x*x[i])
             a += np.asmatrix(xx)
        if ahodo[npt+1] is 'DECO':
             x = [0.,
                  0.,
                  0.,
                  0.,
                  np.cos(aha[npt+1]*dtor),
                  np.sin(aha[npt+1]*dtor),
                  1.,
                - np.cos(lat*dtor)*np.cos(aha[npt+1]*dtor)*np.sin(adec[npt+1]*dtor) 
                + np.sin(lat*dtor)*np.cos(adec[npt+1]*dtor)]
             if nparm == 10:
                 x.append(0.)
                 x.append(adec[npt+1]-lat)
             x = np.array(x)
             n += 1
             b += np.array(apo[npt+1]*x)
             xx = []
             for i in range(len(x)):
                 xx.append(x*x[i])
             a += np.asmatrix(xx)

        npt += 2

    if npt == 0:
        print 'EQ_MOUNTCAL: File read error.'
        return

    # Solve equation for pointing parameters
    p = np.linalg.solve(a,b)

    fit = []
    hafit = []
    decfit = []
    # Do fit to all of the measurements using the solution for P
    #   HAO = P1 - P2 cos(LAT)*cos(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) 
    #            - P6 cos(HA)*tan(DEC) + P9 HA
    #   DECO = P5 cos(HA) + P6 sin(HA) + P7 + P8 [cos(LAT)*cos(HA)*sin(DEC) + sin(LAT)*cos(DEC)]
    #            + P10 DEC
    for i in range(npt):
        if ahodo[i] is 'HAO':
            ha = p[0] - p[1]*np.cos(lat*dtor)*np.sin(aha[i]*dtor)/np.cos(adec[i]*dtor) \
                      + p[2]*np.tan(adec[i]*dtor) \
                      - p[3]/np.cos(adec[i]*dtor) \
                      + p[4]*np.sin(aha[i]*dtor)*np.tan(adec[i]*dtor) \
                      - p[5]*np.cos(aha[i]*dtor)*np.tan(adec[i]*dtor) #+ sfac*aha[i]
            if nparm == 10:
                ha = ha + p[8]*aha[i]
            fit.append(ha)
            hafit.append(ha)
        elif ahodo[i] is 'DECO':
            dec = p[4]*np.cos(aha[i]*dtor) + p[5]*np.sin(aha[i]*dtor) \
                      + p[6] \
                      - p[7]*(np.cos(lat*dtor)*np.cos(aha[i]*dtor)*np.sin(adec[i]*dtor) \
                            - np.sin(lat*dtor)*np.cos(adec[i]*dtor)) #+ sfac*(adec[i]-lat)
            if nparm == 10:
                dec = dec + p[9]*(adec[i]-lat)
            fit.append(dec)
            decfit.append(dec)

    # Calculate the difference between the measured offsets and the fitted ones
    diff = np.array(apo) - np.array(fit)

    # Print results
    print ' '
    print '\     HA    DEC      MEASURED   FITTED  DIFFERENCE (deg)'
    print ' '
    for i in range(npt):
        print "\ {:6.1f}{:6.1f}   {:s}= {:7.3f}  {:7.3f}  {:7.3f}".format(aha[i],adec[i],ahodo[i],apo[i],fit[i],diff[i])

    # Calculate an appropriate residual
    rmsum = diff.std()
    origsum = (apo**2).sum()
    origsum = np.sqrt(origsum/(n-nparm))

    # If these are optical stellar measurements, change the sign of the corrections,
    # since they are positions of the center of the field, not the offset of the source
    # from the center.
    if star: p = -p

    # Print residual and solution
    print ' '
    print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
    print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'
    
    print ''.join('{:7d}'.format(int(int(aligntab[k])+v*10000)) for k,v in enumerate(p)),' ALIGNPARM \ Updated pointing parameters'

#    pha = Hour Angle list for plotting
#    pdec = Declination list for plotting
#    hapo = HA pointing offset for plotting
#    decpo = Dec pointing offset for plotting
#    hafit = Fit for HA for plotting
#    decfit = Fit for Declination for plotting
    # Hour Angle -- get index for sort by ascending order
    print 'start of plot'
    pha = np.array(pha)
    ind = pha.argsort()
    hapo = np.array(hapo)
    decpo = np.array(decpo)
    print 'setting drange'
    drange = np.sqrt((hapo.max() - hapo.min())**2 + (decpo.max() - decpo.min())**2)
    matplotlib.rcParams.update({'font.size':12})
    plt.subplot(221)
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('HA Offset [deg]')
    plt.title('HA Offsets and Fit')
    plt.plot(pha[ind],np.array(hapo)[ind],'o')
    plt.plot(pha[ind],np.array(hafit)[ind])

    plt.subplot(222)
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Declination Offsets and Fit')
    plt.plot(pha[ind],np.array(decpo)[ind],'o')
    plt.plot(pha[ind],np.array(decfit)[ind])

    plt.subplot(223)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-60,60,-25,45])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('Sky Coverage')
    plt.plot(pha[ind],np.array(pdec)[ind],'+')

    plt.subplot(224)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-0.25,0.25,-0.25,0.25])
    plt.xlabel('Hour Angle Offset [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Pointing Relative to Solar Disk')
    th = np.linspace(0,2*np.pi,100)
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),np.array(hapo),np.array(decpo),'+')
    # Set plot size
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    tok = filename.split('/')
    stem = tok[len(tok)-1].split('.')[0]
    plt.savefig(stem+'.pdf')
    plt.close()
    
def anta_mountcal(filename=None, param_string=None, star=False, stepsize=None, sfac=0.0):

    lat = 37.233170       # OVSA Latitude (degrees)
    nparm = 8
    # Open the file containing the raw pointing data
    if filename is None:
        print 'Error: Must specify an input file name.'
        return
    if stepsize != None:
        nparm = 10
    if param_string is None:
        if nparm == 10:
            param_string = '0 0 0 0 0 0 0 0 0 0'
        else:
            param_string = '0 0 0 0 0 0 0 0'

    f = open(filename,'r')

    # Read first line to get old alignment parameters
    # line = f.readline()
    aligntab = param_string.strip().split() # Read old alignment parameters into ALIGNTAB
    # Read rest of file into lines
    inlines = f.readlines()
    f.close()
    lines = []
    for line in inlines:
        if line[0] != '#':
            lines.append(line)
    nlines = len(lines)

##    !p.multi = [0,2,1,0,0]
##    !p.charsize = 1

##    saveresult = ' '

     # Set some initializations
    n = 0          # Number of "good" measurements
    npt = 0        # Total number of entries (some could be "flagged" bad)

    x = np.zeros((nparm,nparm))
    a = np.asmatrix(x)   # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)  # Array to create from measurements
    aha = []
    phah = []
    pdech = []
    phad = []
    pdecd = []
    hapo = []
    decpo = []
    adec = []
    apo = []
    ahodo = []
    hodo = []
    name = []
    dtor = np.pi/180.

    for line in lines:
        nam = line[:8]
        info = map(float,line.strip().split()[3:7])
        ha,dec,dra,ddec = info
        ha = ha#*180/np.pi
        dec = dec#*180/np.pi

        if not np.isnan(dra):
            # Enter line contents into arrays for later use
            ahodo.append('HAO')
            aha.append(ha)
            adec.append(dec)
            phah.append(ha)
            pdech.append(dec)
            # ael[npt] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
            apo.append(-dra)#/np.cos(dec*np.pi/180.)     # Change sign for RA -> HA
            hapo.append(-dra)#/np.cos(dec*np.pi/180.))
            # Coordinate parameters are contained in variable X, which are different
            # for azimuth and elevation measurements
            x = [1.,
                 -np.cos(lat*dtor)*np.sin(aha[npt]*dtor)/np.cos(adec[npt]*dtor),
                 np.tan(adec[npt]*dtor),
                 -1./np.cos(adec[npt]*dtor),
                 np.sin(aha[npt]*dtor)*np.tan(adec[npt]*dtor),
                 -np.cos(aha[npt]*dtor)*np.tan(adec[npt]*dtor),
                 0.,
                 0.]
            if nparm == 10:
                 x.append(aha[npt])
                 x.append(0.)
            x = np.array(x)
            n += 1
            b += np.array(apo[npt]*x)
            xx = []
            for i in range(len(x)):
                xx.append(x*x[i])
            a += np.asmatrix(xx)
            npt += 1

        if not np.isnan(ddec):
            ahodo.append('DECO')
            aha.append(ha)
            adec.append(dec)
            phad.append(ha)
            pdecd.append(dec)
            #ael[npt+1] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
            apo.append(ddec)
            decpo.append(ddec)
        
            x = [0.,
                 0.,
                 0.,
                 0.,
                 np.cos(aha[npt]*dtor),
                 np.sin(aha[npt]*dtor),
                  1.,
                + np.cos(lat*dtor)*np.cos(aha[npt]*dtor)*np.sin(adec[npt]*dtor) 
                - np.sin(lat*dtor)*np.cos(adec[npt]*dtor)]
            if nparm == 10:
                x.append(0.)
                x.append(adec[npt]-lat)
            x = np.array(x)
            n += 1
            b += np.array(apo[npt]*x)
            xx = []
            for i in range(len(x)):
                xx.append(x*x[i])
            a += np.asmatrix(xx)
            npt += 1

    if npt == 0:
        print 'EQ_MOUNTCAL: File read error.'
        return

    # Solve equation for pointing parameters
    p = np.linalg.solve(a,b)

    fit = []
    hafit = []
    decfit = []
    # Do fit to all of the measurements using the solution for P
    #   HAO = P1 - P2 cos(LAT)*cos(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) 
    #            - P6 cos(HA)*tan(DEC) + P9 HA
    #   DECO = P5 cos(HA) + P6 sin(HA) + P7 + P8 [cos(LAT)*cos(HA)*sin(DEC) + sin(LAT)*cos(DEC)]
    #            + P10 DEC
    for i in range(npt):
        if ahodo[i] is 'HAO':
            ha = p[0] - p[1]*np.cos(lat*dtor)*np.sin(aha[i]*dtor)/np.cos(adec[i]*dtor) \
                      + p[2]*np.tan(adec[i]*dtor) \
                      - p[3]/np.cos(adec[i]*dtor) \
                      + p[4]*np.sin(aha[i]*dtor)*np.tan(adec[i]*dtor) \
                      - p[5]*np.cos(aha[i]*dtor)*np.tan(adec[i]*dtor) #+ sfac*aha[i]
            if nparm == 10:
                ha = ha + p[8]*aha[i]
            fit.append(ha)
            hafit.append(ha)
        elif ahodo[i] is 'DECO':
            dec = p[4]*np.cos(aha[i]*dtor) + p[5]*np.sin(aha[i]*dtor) \
                      + p[6] \
                      + p[7]*(np.cos(lat*dtor)*np.cos(aha[i]*dtor)*np.sin(adec[i]*dtor) \
                            - np.sin(lat*dtor)*np.cos(adec[i]*dtor)) #+ sfac*(adec[i]-lat)
            if nparm == 10:
                dec = dec + p[9]*(adec[i]-lat)
            fit.append(dec)
            decfit.append(dec)

    # Calculate the difference between the measured offsets and the fitted ones
    diff = np.array(apo) - np.array(fit)

    # Print results
    print ' '
    print '\     HA    DEC      MEASURED   FITTED  DIFFERENCE (deg)'
    print ' '
    for i in range(npt):
        print "\ {:6.1f}{:6.1f}   {:s}= {:7.3f}  {:7.3f}  {:7.3f}".format(aha[i],adec[i],ahodo[i],apo[i],fit[i],diff[i])

    # Calculate an appropriate residual
    rmsum = diff.std()
    origsum = (np.array(apo)**2).sum()
    origsum = np.sqrt(origsum/(n-nparm))

    # If these are optical stellar measurements, change the sign of the corrections,
    # since they are positions of the center of the field, not the offset of the source
    # from the center.
    if star: p = -p

    # Print residual and solution
    print ' '
    print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
    print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'
    
    print ''.join('{:7d}'.format(int(int(aligntab[k])+v*10000)) for k,v in enumerate(p)),' ALIGNPARM \ Updated pointing parameters'

#    pha = Hour Angle list for plotting
#    pdec = Declination list for plotting
#    hapo = HA pointing offset for plotting
#    decpo = Dec pointing offset for plotting
#    hafit = Fit for HA for plotting
#    decfit = Fit for Declination for plotting
    # Hour Angle -- get index for sort by ascending order
    print 'start of plot'
    pha = np.array(phah)
    ind = pha.argsort()
    hapo = np.array(hapo)
    decpo = np.array(decpo)
    print 'setting drange'
    drange = np.sqrt((hapo.max() - hapo.min())**2 + (decpo.max() - decpo.min())**2)
    matplotlib.rcParams.update({'font.size':12})
    plt.subplot(221)
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('HA Offset [deg]')
    plt.title('HA Offsets and Fit')
    plt.plot(pha[ind],np.array(hapo)[ind],'o')
    plt.plot(pha[ind],np.array(hafit)[ind])

    plt.subplot(222)
    pha = np.array(phad)
    ind = pha.argsort()
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Declination Offsets and Fit')
    plt.plot(pha[ind],np.array(decpo)[ind],'o')
    plt.plot(pha[ind],np.array(decfit)[ind])

    plt.subplot(223)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-60,60,-25,45])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('Sky Coverage')
    plt.plot(phah,np.array(pdech),'b+')
    plt.plot(phad,np.array(pdecd),'b+')

    plt.subplot(224)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-0.25,0.25,-0.25,0.25])
    plt.xlabel('Hour Angle Offset [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Pointing Relative to Solar Disk')
    th = np.linspace(0,2*np.pi,100)
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),np.array(hapo),np.array(decpo),'+')
    plt.plot(np.array(hapo-hafit),np.array(decpo-decfit),'*')
    # Set plot size
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    tok = filename.split('.')
    #stem = tok[len(tok)-1].split('.')[0]
    plt.savefig(tok[0]+'.pdf')
    plt.close()
    
def anta_broken_mountcal(filename=None, param_string=None, star=False, stepsize=None, sfac=0.0):

    from util import common_val_idx

    lat = 37.233170       # OVSA Latitude (degrees)
    nparm = 8
    # Open the file containing the raw pointing data
    if filename is None:
        print 'Error: Must specify an input file name.'
        return
    if stepsize != None:
        nparm = 10
    if param_string is None:
        if nparm == 10:
            param_string = '0 0 0 0 0 0 0 0 0 0'
        else:
            param_string = '0 0 0 0 0 0 0 0'

    f = open(filename,'r')

    # Read first line to get old alignment parameters
    # line = f.readline()
    aligntab = param_string.strip().split() # Read old alignment parameters into ALIGNTAB
    # Read rest of file into lines
    inlines = f.readlines()
    f.close()
    lines = []
    for line in inlines:
        if line[0] != '#':
            lines.append(line)
    nlines = len(lines)

##    !p.multi = [0,2,1,0,0]
##    !p.charsize = 1

##    saveresult = ' '

     # Set some initializations
    n = 0          # Number of "good" measurements
    npt = 0        # Total number of entries (some could be "flagged" bad)

    x = np.zeros((nparm,nparm))
    a = np.asmatrix(x)   # Matrix to create from coordinates   A#B = P
    b = np.zeros(nparm)  # Array to create from measurements
    aha = []
    phah = []
    pdech = []
    phad = []
    pdecd = []
    hapo = []
    decpo = []
    adec = []
    apo = []
    ahodo = []
    hodo = []
    name = []
    hidx = []     # Index of non-nan lines for HAO
    didx = []     # Index of non-nan lines for DECO
    dtor = np.pi/180.

    for k,line in enumerate(lines):
        nam = line[:8]
        info = map(float,line.strip().split()[3:7])
        ha,dec,dra,ddec = info
        ha = ha#*180/np.pi
        dec = dec#*180/np.pi

        if not np.isnan(dra):
            # Enter line contents into arrays for later use
            ahodo.append('HAO')
            aha.append(ha)
            adec.append(dec)
            phah.append(ha)
            pdech.append(dec)
            hidx.append(k)
            # ael[npt] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
            apo.append(-dra)#/np.cos(dec*np.pi/180.)     # Change sign for RA -> HA
            hapo.append(-dra)#/np.cos(dec*np.pi/180.))
            # Coordinate parameters are contained in variable X, which are different
            # for azimuth and elevation measurements
            x = [1.,
                 -np.cos(lat*dtor)*np.sin(aha[npt]*dtor)/np.cos(adec[npt]*dtor),
                 np.tan(aha[npt]*dtor),
                 -1./np.cos(adec[npt]*dtor),
                 np.sin(aha[npt]*dtor)*np.tan(aha[npt]*dtor),
                 -np.cos(aha[npt]*dtor)*np.tan(aha[npt]*dtor),
                 0.,
                 0.]
            if nparm == 10:
                 x.append(aha[npt])
                 x.append(0.)
            x = np.array(x)
            n += 1
            b += np.array(apo[npt]*x)
            xx = []
            for i in range(len(x)):
                xx.append(x*x[i])
            a += np.asmatrix(xx)
            npt += 1

        if not np.isnan(ddec):
            ahodo.append('DECO')
            aha.append(ha)
            adec.append(dec)
            phad.append(ha)
            pdecd.append(dec)
            #ael[npt+1] = el+(1/60.)*(0.0019279 + 1.02/np.tan((el + 10.3/(el+5.1))*dtor))  # refraction corr.
            apo.append(ddec)
            decpo.append(ddec)
            didx.append(k)
        
            x = [0.,
                 0.,
                 0.,
                 0.,
                 np.cos(aha[npt]*dtor),
                 np.sin(aha[npt]*dtor),
                  1.,
                + np.cos(lat*dtor)*np.cos(aha[npt]*dtor)*np.sin(adec[npt]*dtor) 
                - np.sin(lat*dtor)*np.cos(adec[npt]*dtor)]
            if nparm == 10:
                x.append(0.)
                x.append(adec[npt]-lat)
            x = np.array(x)
            n += 1
            b += np.array(apo[npt]*x)
            xx = []
            for i in range(len(x)):
                xx.append(x*x[i])
            a += np.asmatrix(xx)
            npt += 1

    if npt == 0:
        print 'EQ_MOUNTCAL: File read error.'
        return

    # Solve equation for pointing parameters
    p = np.linalg.solve(a,b)

    fit = []
    hafit = []
    decfit = []
    # Do fit to all of the measurements using the solution for P
    #   HAO = P1 - P2 cos(LAT)*cos(HA)*sec(DEC) + P3 tan(DEC) - P4 sec(DEC) + P5 sin(HA)*tan(DEC) 
    #            - P6 cos(HA)*tan(DEC) + P9 HA
    #   DECO = P5 cos(HA) + P6 sin(HA) + P7 + P8 [cos(LAT)*cos(HA)*sin(DEC) + sin(LAT)*cos(DEC)]
    #            + P10 DEC
    for i in range(npt):
        if ahodo[i] is 'HAO':
            ha = p[0] - p[1]*np.cos(lat*dtor)*np.sin(aha[i]*dtor)/np.cos(adec[i]*dtor) \
                      + p[2]*np.tan(aha[i]*dtor) \
                      - p[3]/np.cos(adec[i]*dtor) \
                      + p[4]*np.sin(aha[i]*dtor)*np.tan(aha[i]*dtor) \
                      - p[5]*np.cos(aha[i]*dtor)*np.tan(aha[i]*dtor) #+ sfac*aha[i]
            if nparm == 10:
                ha = ha + p[8]*aha[i]
            fit.append(ha)
            hafit.append(ha)
        elif ahodo[i] is 'DECO':
            dec = p[4]*np.cos(aha[i]*dtor) + p[5]*np.sin(aha[i]*dtor) \
                      + p[6] \
                      + p[7]*(np.cos(lat*dtor)*np.cos(aha[i]*dtor)*np.sin(adec[i]*dtor) \
                            - np.sin(lat*dtor)*np.cos(adec[i]*dtor)) #+ sfac*(adec[i]-lat)
            if nparm == 10:
                dec = dec + p[9]*(adec[i]-lat)
            fit.append(dec)
            decfit.append(dec)

    # Calculate the difference between the measured offsets and the fitted ones
    diff = np.array(apo) - np.array(fit)

    # Print results
    print ' '
    print '\     HA    DEC      MEASURED   FITTED  DIFFERENCE (deg)'
    print ' '
    for i in range(npt):
        print "\ {:6.1f}{:6.1f}   {:s}= {:7.3f}  {:7.3f}  {:7.3f}".format(aha[i],adec[i],ahodo[i],apo[i],fit[i],diff[i])

    # Calculate an appropriate residual
    rmsum = diff.std()
    origsum = (np.array(apo)**2).sum()
    origsum = np.sqrt(origsum/(n-nparm))

    # If these are optical stellar measurements, change the sign of the corrections,
    # since they are positions of the center of the field, not the offset of the source
    # from the center.
    if star: p = -p

    # Print residual and solution
    print ' '
    print '\ RMS Residual={:5.3f}   RMS Original={:5.3f}'.format(rmsum,origsum)
    print '\\',''.join('{:8.4f}'.format(k) for k in p), '<- UPDATE'
    
    print ''.join('{:7d}'.format(int(int(aligntab[k])+v*10000)) for k,v in enumerate(p)),' ALIGNPARM \ Updated pointing parameters'

#    pha = Hour Angle list for plotting
#    pdec = Declination list for plotting
#    hapo = HA pointing offset for plotting
#    decpo = Dec pointing offset for plotting
#    hafit = Fit for HA for plotting
#    decfit = Fit for Declination for plotting
    # Hour Angle -- get index for sort by ascending order
    print 'start of plot'
    pha = np.array(phah)
    ind = pha.argsort()
    hapo = np.array(hapo)
    decpo = np.array(decpo)
    print 'setting drange'
    drange = np.sqrt((hapo.max() - hapo.min())**2 + (decpo.max() - decpo.min())**2)
    matplotlib.rcParams.update({'font.size':12})
    plt.subplot(221)
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('HA Offset [deg]')
    plt.title('HA Offsets and Fit')
    plt.plot(pha[ind],np.array(hapo)[ind],'o')
    plt.plot(pha[ind],np.array(hafit)[ind])

    plt.subplot(222)
    pha = np.array(phad)
    ind = pha.argsort()
    if drange > 0.25:
        plt.axis([-60,60,-3,3])
    else:
        plt.axis([-60,60,-0.25,0.25])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Declination Offsets and Fit')
    plt.plot(pha[ind],np.array(decpo)[ind],'o')
    plt.plot(pha[ind],np.array(decfit)[ind])

    plt.subplot(223)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-60,60,-25,45])
    plt.xlabel('Hour Angle [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('Sky Coverage')
    plt.plot(phah,np.array(pdech),'b+')
    plt.plot(phad,np.array(pdecd),'b+')

    # Can only plot 2D pointing for pairs of points, so use hidx and didx to find them:
    hgood, dgood = common_val_idx(hidx, didx)
    plt.subplot(224)
    plt.axis('equal')
    plt.axis('scaled')
    plt.axis([-0.25,0.25,-0.25,0.25])
    plt.xlabel('Hour Angle Offset [deg]')
    plt.ylabel('Declination Offset [deg]')
    plt.title('Pointing Relative to Solar Disk')
    th = np.linspace(0,2*np.pi,100)
    plt.plot(0.25*np.cos(th),0.25*np.sin(th),color='orange')
    for i in range(len(hgood)):
        plt.plot(np.array(hapo[hgood[i]]),np.array(decpo[dgood[i]]),'b+')
        plt.plot(np.array(hapo-hafit)[hgood[i]],np.array(decpo-decfit)[dgood[i]],'g*')
    # Set plot size
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    tok = filename.split('.')
    #stem = tok[len(tok)-1].split('.')[0]
    plt.savefig(tok[0]+'.pdf')
    plt.close()
    
def ptg_model(ha, dec, param_string='-1875 1442 1470 -162 862 -127 9878 468',alt=False):
    ''' Given an HA and Dec [degrees], and a string of pointing coefficients, 
        returns the HA and Dec offsets to be added to the current HA and Dec.
    '''
    lat = 37.233170
    dtor = np.pi/180.
    p = np.array(map(float,param_string.strip().split()))/10000. # Read old alignment parameters into ALIGNTAB

    x = np.array([1.,
           -np.cos(lat*dtor)*np.sin(ha*dtor)/np.cos(dec*dtor),
            np.tan(dec*dtor),
           -1./np.cos(dec*dtor),
            np.sin(ha*dtor)*np.tan(dec*dtor),
           -np.cos(ha*dtor)*np.tan(dec*dtor),
            0.,
            0.])

    if alt:
        x = np.array([1.,
           -np.cos(lat*dtor)*np.sin(ha*dtor)/np.cos(dec*dtor),
            np.tan(ha*dtor),
           -1./np.cos(dec*dtor),
            np.sin(ha*dtor)*np.tan(ha*dtor),
           -np.cos(ha*dtor)*np.tan(ha*dtor),
            0.,
            0.])

    y = np.array([0.,
               0.,
               0.,
               0.,
               np.cos(ha*dtor),
               np.sin(ha*dtor),
               1.,
               np.cos(lat*dtor)*np.cos(ha*dtor)*np.sin(dec*dtor) 
               - np.sin(lat*dtor)*np.cos(dec*dtor)])
         
    return sum(p*x), sum(p*y)
