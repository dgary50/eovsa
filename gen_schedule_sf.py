# History:
#   2016-Mar-08  DG
#      Started this history log.  Made maximum delay 16000, to reflect
#      new maximum in the 16-ant correlator design.
#   2017-Mar-05  DG
#      Made maximum delay 32000, to reflect new maximum in the now-working
#      300 MHz correlator.
#
import struct
import numpy as np

def gen_schedule_sf(sf_dict,mk_xml=False):
    '''Writes the schedule stateframe items from the stateframe dictionary 
       created by the schedule. Optionally creates the corresponding XML
       file but always returns its file name (in xmlfile in the /tmp directory).
       Also returns the binary data buffer (buf), and format string (fmt).
       The format string fmt can be used with struct.unpack() to read 
       the data from buf, although that usage is not anticipated except 
       for testing.

       This routine does something sensible even if the supplied
       dictionary sf_dict is empty (i.e. is {}).
    '''
    dtor = np.pi/180.
    xmlfile = '/tmp/schedule_stateframe.xml'
    if mk_xml:
        xml = open(xmlfile,'w')
        xml.write('<Cluster>\n')
        xml.write('<Name>Data</Name>\n')
        xml.write('<NumElts>13</NumElts>\n')

    # Schedule_Timestamp (double) [s, in LabVIEW format]
    # To be compatible with other timestamps in the stateframe, this
    # will be in LabVIEW format, which is s since 1904/01/01 (don't ask).
    # It is the time (should be exact second, no microseconds) for
    # which the UVW coordinates and Delays are calculated.
    item = sf_dict.get('timestamp',0.0)
    fmt = '<d'
    buf = struct.pack('d',item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>Timestamp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')
    
    # Schedule version (double) [N/A]
    # Version of the schedule stateframe.
    item = sf_dict.get('Version',0.4)
    fmt += 'd'
    buf += struct.pack('d',item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>Version</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')
    
    # Scan_State (unsigned integer bool)
    # Flag (=1 to indicate that DPP should be recording data, =0 otherwise)
    item = sf_dict.get('scan_state',0)
    fmt += 'i'
    buf += struct.pack('i',item)
    if mk_xml:
        xml.write('<I32>\n')
        xml.write('<Name>ScanState</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</I32>\n')

    # Phase_Tracking (unsigned integer bool)
    # Flag (=1 to indicate that uvw coordinates are valid, =0 otherwise)
    item = sf_dict.get('phase_tracking',0)
    fmt += 'I'
    buf += struct.pack('I',item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>PhaseTracking</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # UVW (3 x 16 array of doubles) [ns]
    # u, v, w coordinates for each antenna, relative to antenna 1.
    # Default is array of zeros (=> not tracking phase center)
    item = sf_dict.get('uvw',np.array([[0.0,0.0,0.0]]*16))
    # Write dimensions into data stream
    fmt += 'II'
    buf += struct.pack('II',3,16)
    fmt += str(3*16)+'d'
    for i in range(16):
        buf += struct.pack('3d',item[i,0],item[i,1],item[i,2])
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>UVW</Name>\n')
       xml.write('<Dimsize>3</Dimsize><Dimsize>16</Dimsize>\n<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
       xml.write('</Array>\n')

    # Delay (length 16 x 2 array of doubles) [ns]
    # Geometric delay (-w coordinate) for each antenna, relative to antenna 1,
    # for current time (stateframe timestamp), and again for current time plus
    # 1 s (delay1).
    # Default is array of zeros (=> not tracking phase center)
    # Write dimensions into data stream
    fmt += 'II'
    buf += struct.pack('II',16,2)
    item = sf_dict.get('delay',np.zeros(16))
    fmt += '32d'
    for i in item:
        buf += struct.pack('d',i)
    item = sf_dict.get('delay1',np.zeros(16))
    for i in item:
        buf += struct.pack('d',i)
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>Delay</Name>\n')
       xml.write('<Dimsize>16</Dimsize><Dimsize>2</Dimsize>\n<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
       xml.write('</Array>\n')

    # Azimuth (length 16 array of doubles) [radians]
    # Current azimuth pointing for each antenna
    # Default is zero for any missing antennas
    fmt += 'I'
    buf += struct.pack('I',16)
    item = sf_dict.get('ActualAzimuth',np.zeros(15))*dtor  # Convert degrees to radians
    if len(item) == 15:
        item = np.append(item,0.0)  # Add antenna 16, whose Az position is zero by definition
    fmt += '16d'
    for i in item:
        buf += struct.pack('d',i)
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>Azimuth</Name>\n')
       xml.write('<Dimsize>16</Dimsize>\n<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
       xml.write('</Array>\n')

    # Elevation (length 16 array of doubles) [radians]
    # Current elevation pointing for each antenna
    # Default is zero for any missing antennas
    fmt += 'I'
    buf += struct.pack('I',16)
    item = sf_dict.get('ActualElevation',np.zeros(15))*dtor  # Convert degrees to radians
    if len(item) == 15:
        item = np.append(item,0.0)  # Add antenna 16, whose El position is zero by definition
    fmt += '16d'
    for i in item:
        buf += struct.pack('d',i)
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>Elevation</Name>\n')
       xml.write('<Dimsize>16</Dimsize>\n<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
       xml.write('</Array>\n')

    # Chi (length 16 array of doubles) [radians]
    # Current parallactic angle for each antenna
    # Default is zero for any missing antennas
    fmt += 'I'
    buf += struct.pack('I',16)
    item = sf_dict.get('ParallacticAngle',np.zeros(15))*dtor  # Convert degrees to radians
    item = np.append(item,0.0)  # Add antenna 16, whose parallactic angle is zero by definition
    fmt += '16d'
    for i in item:
        buf += struct.pack('d',i)
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>Chi</Name>\n')
       xml.write('<Dimsize>16</Dimsize>\n<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
       xml.write('</Array>\n')
 
    # Antenna Tracking (unsigned integer bool array)
    # Flags (=1 to indicate that antenna is tracking)
    item = sf_dict.get('TrackFlag',np.array([False]*16))
    # Write dimension (16) into data stream
    fmt += 'I'
    buf += struct.pack('I',16)
    fmt += '16I'
    for i in item:
        buf += struct.pack('I',i)
    if mk_xml:
       xml.write('<Array>\n')
       xml.write('<Name>TrackFlag</Name>\n')
       xml.write('<Dimsize>16</Dimsize>\n<U32>\n<Name></Name>\n<Val></Val>\n</U32>\n')
       xml.write('</Array>\n')

    # XML encapsulation of weather info
    if mk_xml:
       xml.write('<Cluster>\n')
       xml.write('<Name>Weather</Name>\n')
       xml.write('<NumElts>10</NumElts>\n')

    # Wind (float) [mph]
    # Current "instantaneous" wind speed
    # Default is zero
    # item = sf_dict.get('Wind Speed',0.0)
    item = sf_dict.get('mtWindSpeed',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>Wind</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Wind_Direction (float) [radians]
    # Current "instantaneous" wind direction +ve East of North
    # Default is zero
    # item = sf_dict.get('Adjusted Wind Direction',0.0)
    item = sf_dict.get('mtAdjWindDir',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    item = item*dtor
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>WindDirection</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # AvgWind (float) [mph]
    # 10-min rolling average wind speed
    # Default is zero
    # item = sf_dict.get('2m Wind Speed',0.0)
    item = sf_dict.get('mt2MinRollAvgWindSpeed',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>AvgWind</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # AvgWind_Direction (float) [radians]
    # 10-min rolling average wind direction +ve East of North
    # Default is zero
    # item = sf_dict.get('2m Wind Direction',0.0)
    item = sf_dict.get('mt2MinRollAvgWindDir',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    item = item*dtor
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>AvgWindDirection</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # AvgWind_Gust (float) [mph]
    # 10-min rolling max wind speed
    # Default is zero
    # item = sf_dict.get('10m Wind Gust Speed',0.0)
    item = sf_dict.get('mt10MinWindGustSpeed',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    item = item
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>AvgWindGust</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Temperature (float) [degrees F]
    # Current air temperature
    # Default is zero
    # item = sf_dict.get('Temperature',0.0)
    item = sf_dict.get('mtTemp1',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>Temperature</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Pressure (float) [mbar]
    # Current air pressure
    # Default is zero
    # item = sf_dict.get('Raw Barometric Pressure',0.0)
    item = sf_dict.get('mtRawBaromPress',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>Pressure</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Humidity (float) [percent]
    # Current relative humidity
    # Default is zero
    # item = sf_dict.get('Relative Humidity',0.0)
    item = sf_dict.get('mtRelHumidity',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>Humidity</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Rain_Rate (float) [inches/hour]
    # Running 5-min rate of rainfall
    # Default is zero
    # item = sf_dict.get('Rain Rate',0.0)
    item = sf_dict.get('mtRainRate',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>RainRate</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    # Rain_Today (float) [inches]
    # Rainfall in past 24-h
    # Default is zero
    # item = sf_dict.get('Rain Today',0.0)
    item = sf_dict.get('mtRainToday',0.0)
    try:
        item = float(item)
    except ValueError:
        item = 0.0
    fmt += 'f'
    buf += struct.pack('f',item)
    if mk_xml:
        xml.write('<SGL>\n')
        xml.write('<Name>RainToday</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')

    if mk_xml:
        # End Weather Cluster
        xml.write('</Cluster>\n')

    # XML encapsulation of solar power station info (array of 2)
    if mk_xml:
        xml.write('<Array>\n')
        xml.write('<Name>SolarPower</Name>\n')
        xml.write('<Dimsize>2</Dimsize>\n')
    # Put array dimension into data
    fmt += 'I'
    buf += struct.pack('I',2)
    # This is an array, but the definitions only go in once
    # Here are the definitions for the 8 elements in the cluster
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name/>\n')
        xml.write('<NumElts>8</NumElts>\n')
        xml.write('<DBL>\n')
        xml.write('<Name>Timestamp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Charge</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Volts</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Amps</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<I32>\n')
        xml.write('<Name>AmpHours</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</I32>\n')
        xml.write('<I32>\n')
        xml.write('<Name>BatteryTemp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</I32>\n')
        xml.write('<I32>\n')
        xml.write('<Name>TransformerTemp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</I32>\n')
        xml.write('<I32>\n')
        xml.write('<Name>FETTemp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</I32>\n')
        # End SolarPower Cluster
        xml.write('</Cluster>\n')

    # Get the solar power array of dictionaries (defaults to two empty dicts)
    solpwr = sf_dict.get('SolPwr',[{},{}])
    # The actual data has to go in twice
    for i in range(2):
        # Solar Power Timestamp (double) [s, in LabVIEW format]
        # Time that data were reported, in seconds since 1904 Jan 1
        # Default is zero
        item = solpwr[i].get('Time',0.0)
        fmt += 'd'
        buf += struct.pack('d',item)

        # Charge (unsigned int) [% of full charge]
        # Default is zero
        item = solpwr[i].get('Charge',0)
        fmt += 'I'
        buf += struct.pack('I',item)

        # Volts (float) [DC voltage]
        # Default is zero
        item = solpwr[i].get('Volts',0.0)
        fmt += 'f'
        buf += struct.pack('f',item)

        # Amps (float) [DC amps]
        # Default is zero
        item = solpwr[i].get('Amps',0.0)
        fmt += 'f'
        buf += struct.pack('f',item)

        # AmpHours (int) [Charge available]
        # Default is zero
        item = solpwr[i].get('AmpHours',0.0)
        fmt += 'i'
        buf += struct.pack('i',item)

        # BatteryTemp (int) [degrees C]
        # Battery Temperature in degrees C
        # Default is zero
        item = solpwr[i].get('BatteryTemp',0)
        fmt += 'i'
        buf += struct.pack('i',item)

        # TransformerTemp (int) [degrees C]
        # Transformer Temperature in degrees C
        # Default is zero
        item = solpwr[i].get('TransformerTemp',0)
        fmt += 'i'
        buf += struct.pack('i',item)

        # FETTemp (int) [degrees C]
        # FET Temperature in degrees C
        # Default is zero
        item = solpwr[i].get('FETTemp',0)
        fmt += 'i'
        buf += struct.pack('i',item)

    if mk_xml:
        # End SolarPower Array
        xml.write('</Array>\n')

    # XML encapsulation of roach info (array of 8)
    if mk_xml:
        xml.write('<Array>\n')
        xml.write('<Name>Roach</Name>\n')
        xml.write('<Dimsize>8</Dimsize>\n')
    # Put array dimension into data
    fmt += 'I'
    buf += struct.pack('I',8)
    # This is an array, but the definitions only go in once
    # Here are the definitions for the 8 elements in the cluster
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name/>\n')
        xml.write('<NumElts>30</NumElts>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Status</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.12v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.1v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.1v5</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.1v8</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.2v5</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.3v3</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Current.5v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Fan.chs0</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Fan.csh1</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Fan.csh2</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Fan.fpga</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Temp.ambient</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Temp.fpga</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Temp.inlet</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Temp.outlet</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Temp.ppc</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.12v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.1v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.1v5</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.1v8</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.2v5</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.3v3</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.3v3aux</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.5v</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<SGL>\n')
        xml.write('<Name>Voltage.5vaux</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</SGL>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Delay0x</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Delay0y</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Delay1x</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        xml.write('<U32>\n')
        xml.write('<Name>Delay1y</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')
        # End Roaches Cluster
        xml.write('</Cluster>\n')

    # Get the roach array of dictionaries (defaults to eight empty dicts)
    sensors = sf_dict.get('sensors',[{},{},{},{},{},{},{},{}])
    delays = sf_dict.get('delays',[{},{},{},{},{},{},{},{}])
    # The actual data has to go in eight times
    for i in range(8):
        # Handle case of empty dictionary
        if sensors[i] == {}:
            fmt += 'IfffffffIIIIffffffffffffff'
            buf += struct.pack('I',0)
            for j in range(7):
                buf += struct.pack('f',0.0)
            for j in range(4):
                buf += struct.pack('I',0)
            for j in range(14):
                buf += struct.pack('f',0.0)
        
        else:
            # Status (int) [bit list]
            # Each bit signifies whether corresponding sensor is nominal [0] or in error [1]
            # Order is alphabetical, msb to lsb.
            # Default is zero
            status = 0
            keys = sorted(sensors[i].keys())
            for key in keys:
                if key.find('status') != -1:
                    if sensors[i][key] == 'nominal':
                        status = status<<1
                    else:
                        status = (status<<1) + 1
            item = status
            fmt += 'I'
            buf += struct.pack('I',item)

            # Just add the sorted list of sensors
            # Default is zero
            for key in keys:
                if key.find('status') == -1:
                    item = sensors[i][key]
                    if key.find('fan') == -1:
                        # This is not a fan, so type is float (SGL)
                        fmt += 'f'
                        buf += struct.pack('f',item)
                    else:
                        # This is a fan, so type is uint (U32)
                        fmt += 'I'
                        buf += struct.pack('I',item)

        # Handle case of empty dictionary
        if delays[i] == {}:
            fmt += 'IIII'
            for j in range(4):
                buf += struct.pack('I',0)
        else:
            # Not empty, so add delays
            # Default is zero
            item = delays[i].get('dx0',0.0) 
            fmt += 'I'
            item = np.clip(item,0,32000)
            buf += struct.pack('I',item)
            item = delays[i].get('dy0',0.0)
            fmt += 'I'
            item = np.clip(item,0,32000)
            buf += struct.pack('I',item)
            item = delays[i].get('dx1',0.0)
            fmt += 'I'
            item = np.clip(item,0,32000)
            buf += struct.pack('I',item)
            item = delays[i].get('dy1',0.0)
            fmt += 'I'
            item = np.clip(item,0,32000)
            buf += struct.pack('I',item)
        
    if mk_xml:
        # End Roach Array
        xml.write('</Array>\n')

    if mk_xml:
        # End Schedule Cluster
        xml.write('</Cluster>\n')
        xml.close()

    return fmt, buf, xmlfile

