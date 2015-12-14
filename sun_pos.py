from math import sin, cos, pi, atan2, atan, asin, tan
def  sun_pos(dd):
    '''This routine is a truncated version of Newcomb's Sun and
       is designed to give apparent angular coordinates (T.E.D) to a
       precision of one second of time.

       Translated from SSW (IDL) routine of the same name
    '''
    dtor = pi/180.
    
    #  Form time in Julian centuries from 1900.0
    t = dd/36525.0

    #  Form sun's mean longitude
    l = (279.696678+((36000.768925*t) % 360.0))*3600.0

    # Allow for ellipticity of the orbit (equation of centre)
    # using the Earth's mean anomoly ME
    me = 358.475844 + ((35999.049750*t) % 360.0)
    ellcor  = (6910.1 - 17.2*t)*sin(me*dtor) + 72.3*sin(2.0*me*dtor)
    l = l + ellcor

    # Allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + ((58517.803875*t) % 360.0) 
    vencorr = 4.8 * cos((299.1017 + mv - me)*dtor) + \
              5.5 * cos((148.3133 +  2.0 * mv  -  2.0 * me )*dtor) + \
              2.5 * cos((315.9433 +  2.0 * mv  -  3.0 * me )*dtor) + \
              1.6 * cos((345.2533 +  3.0 * mv  -  4.0 * me )*dtor) + \
              1.0 * cos((318.15   +  3.0 * mv  -  5.0 * me )*dtor)
    l = l + vencorr

    # Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.529425  +  (( 19139.858500 * t)  %  360.0 )
    marscorr = 2.0 * cos((343.8883 -  2.0 * mm  +  2.0 * me)*dtor ) + \
               1.8 * cos((200.4017 -  2.0 * mm  + me)*dtor)
    l = l + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of
    # Jupiter MJ
    mj = 225.328328  +  (( 3034.6920239 * t)  %  360.0 )
    jupcorr = 7.2 * cos(( 179.5317 - mj + me )*dtor) + \
              2.6 * cos((263.2167  -  mj )*dtor) + \
              2.7 * cos(( 87.1450  -  2.0 * mj  +  2.0 * me)*dtor) + \
              1.6 * cos((109.4933  -  2.0 * mj  +  me )*dtor)
    l = l + jupcorr
    
    # Allow for the Moon's perturbations using the mean elongation of
    # the Moon from the Sun D
    d = 350.7376814  + (( 445267.11422 * t)  %  360.0 )
    mooncorr  = 6.5 * sin(d*dtor)
    l = l + mooncorr

    # Allow for long period terms
    longterm  = + 6.4 * sin(( 231.19  +  20.20 * t )*dtor)
    l  =    l + longterm
    l  =  ( l + 2592000.0)  %  1296000.0 
    longmed = l/3600.0

    # Allow for Aberration
    l  =  l - 20.5

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.183275 - (( 1934.142008 * t ) % 360.0 )
    l  =  l - 17.2 * sin(omega*dtor)

    # Form the True Obliquity
    oblt  = 23.452294 - 0.0130125*t + (9.2*cos(omega*dtor))/3600.0

    # Form Right Ascension and Declination
    l = l/3600.0
    ra  = atan2( sin(l*dtor) * cos(oblt*dtor) , cos(l*dtor) ) / dtor

    if ra < 0.0:
        ra += 360.0

    dec = asin(sin(l*dtor) * sin(oblt*dtor)) / dtor

    return longmed, ra, dec, l, oblt

def  get_pb0r(mjd,arcsec=False):
    '''Given a modified Julian date, return the solar P-angle (degrees),
       B0-angle (degrees), and solar radius (arcmin, or if arcsec=True,
       return solar radius in arcsec)

       Translated from SSW (IDL) routine pb0r().
    '''
    dtor = pi/180.
    de = mjd-15019.5  # Parameters defined starting at noon on 1899/12/31.
#;---------------------------------------------------------------------------
#;  get the longitude of the sun etc.
#;---------------------------------------------------------------------------
    longmed, ra, dec, appl, oblt = sun_pos(de)

#;---------------------------------------------------------------------------
#;  form aberrated longitude
#;---------------------------------------------------------------------------
    lmbda = longmed - (20.5/3600.0)

#;---------------------------------------------------------------------------
#;  form longitude of ascending node of sun's equator on ecliptic
#;---------------------------------------------------------------------------
    node = 73.666666 + (50.25/3600.0)*( (de/365.25) + 50.0 )
    arg = lmbda - node

#;---------------------------------------------------------------------------
#;  calculate P, the position angle of the pole
#;---------------------------------------------------------------------------
    p = (atan(-tan(oblt*dtor) * cos(appl*dtor)) + 
         atan( -0.12722 * cos(arg*dtor))) / dtor

#;---------------------------------------------------------------------------
#;  ... and B0 the tilt of the axis
#;---------------------------------------------------------------------------
    b = asin( 0.12620 * sin(arg*dtor) ) / dtor

#;---------------------------------------------------------------------------
#;  ... and the semi-diameter
#;
#;
#;  Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
#;  and the mean elongation of the Moon from the Sun(D).
#;
#;---------------------------------------------------------------------------
    t = de/36525.0

    mv = 212.6   + ( (58517.80   * t) % 360.0 )
    me = 358.476 + ( (35999.0498 * t) % 360.0 )
    mm = 319.5   + ( (19139.86   * t) % 360.0 )
    mj = 225.3   + ( ( 3034.69   * t) % 360.0 )
    d  = 350.7   + ( (445267.11  * t) % 360.0 )

#;---------------------------------------------------------------------------
#;  Form the geocentric distance(r) and semi-diameter(sd)
#;---------------------------------------------------------------------------
    r = 1.000141 - (0.016748 - 0.0000418*t)*cos(me*dtor) \
      - 0.000140 * cos(2.0*me*dtor)                      \
      + 0.000016 * cos((58.3 + 2.0*mv - 2.0*me)*dtor)    \
      + 0.000005 * cos((209.1 + mv - me)*dtor)           \
      + 0.000005 * cos((253.8 - 2.0*mm + 2.0*me)*dtor)   \
      + 0.000016 * cos(( 89.5 - mj + me)*dtor)           \
      + 0.000009 * cos((357.1 - 2.0*mj + 2.0*me)*dtor)   \
      + 0.000031 * cos(d*dtor)

    sd = (0.2665685/r)*60.0

    if arcsec:
        return p, b, sd*60.
    return p, b, sd

