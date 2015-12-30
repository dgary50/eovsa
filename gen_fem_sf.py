"""
    STARBURST Front End Item Struct Decomposition
    (Based on gen_schedule_sf.py)
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""
# Change-log:
#   2015-12-30  DG
#     Changed RxSel to SelectedRx.  Also swapped FirstStageTemp and
#     SecondStageTemp.  Changed stateframe version accordingly.
#
import numpy as np
import struct
import shutil

# NUMBER OF ELEMENTS IN CLUSTERS:
NELEMENTS_ANTENNA = 6

NELEMENTS_ANTENNA_POWERSTRIP = 10

NELEMENTS_ANTENNA_THERMAL = 9

NELEMENTS_ANTENNA_RECEIVER = 4
NELEMENTS_ANTENNA_RECEIVER_LNA = 6

NELEMENTS_ANTENNA_SERVO = 5
NELEMENTS_ANTENNA_SERVO_AXIS = 7

# POWERSTRIP DEFINITIONS
POWERSTRIP_DEF = ['RFSwitchStatus',
                  'OpticalTxRabbitStatus',
                  'DeltaTauBrickStatus',
                  'ComputerStatus',
                  'LNA12VBiasStatus',
                  'LNA5VBiasBBStatus',
                  'SecondAmpsStatus',
                  'NoiseDiodeStatus']

# THERMAL DEFINITIONS
THERMAL_DEF = ['FirstStageTemp',
               'HiFreq15KPlateTemp',
               'HiFreqLNATemp',
               'HiFreqFeedhornTemp',
               'SecondStageTemp',
               'LowFreqLNATemp',
               'RadiationShieldTemp',
               'LowFreqFeedhornTemp']

# RECEIVER DEFINITIONS
RECEIVER_LNA_DEF = ['DRAINVOLTAGE',
                    'GATEAVOLTAGE',
                    'GATEBVOLTAGE',
                    'DRAINCURRENT',
                    'GATEACURRENT',
                    'GATEBCURRENT']

# AXIS DEFINITIONS
AXIS_DEF = {1: 'ZFocus',
            3: 'PositionAngle',
            4: 'RxSelect'}

# Version Number for FEM stateframe
VERSION = 1.3              # Version Date: 12/30/15
VERSION_DATE = '12.30.15'   # Most recent update (used to write backup file)


def gen_fem_sf(sf_dict, mk_xml=False):
    # Set up file name, format string, and buffer.
    xmlFile = r'tmp/femab_stateframe.xml'
    fmt = '<'
    buf = ''
    xml = None

    # Append XML for antenna clusters. (Note there are two antennas.)
    if mk_xml:
        xml = open(xmlFile, "w")
        xml.write('<Cluster>\n')
        xml.write('<Name>FEM</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA)
                  + '</NumElts>\n')

    append_fmt, append_buf = __antenna(sf_dict, xml, mk_xml)
    fmt += append_fmt
    buf += append_buf

    # Append for end of data cluster
    if mk_xml:
        xml.write('</Cluster>\n')
        xml.close()

        # Make backup copy of XML file
        backup_file = ('starburst/fem_stateframe_v' +
                       str(VERSION) + '_' + VERSION_DATE + '.xml')
        shutil.copyfile(xmlFile, backup_file)

        # Print size of buf
        print 'fem size =', len(buf)
        print 'Modify acc.ini to reflect this if this is a change in size'

    return fmt, buf, xmlFile


def __powerstrip(dict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - PowerStrip:
    # ----------------------------------------------------------------------
    default_statuses = np.zeros(8)
    default_volt = np.zeros(2)
    default_current = np.zeros(2)

    # ----------------------------------------------------------------------
    # XML Cluster setup.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name>PowerStrip</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_POWERSTRIP)
                  + '</NumElts>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 1-8> Status of each device: 0 = off, 1 = on (unsigned int)
    # ----------------------------------------------------------------------
    # Devices are as follows:
    #   0: RF Switch
    #   1: Optical Tx/Rabbit
    #   2: Delta Tau Brick
    #   3: Computer
    #   4: 12V LNA Bias
    #   5: 5V LNA Bias/Beaglebone
    #   6: 2nd Amps (5V)
    #   7: Noise Diode
    # ----------------------------------------------------------------------

    # Pack each status as unsigned int.
    item = dict.get('STATUS', default_statuses)
    for i in range(0, 8):
        fmt += 'I'
        buf += struct.pack('I', int(item[i]))

        # Append to XML file
        if mk_xml:
            xml.write('<U32>\n')
            xml.write('<Name>' + POWERSTRIP_DEF[i] + '</Name>\n')
            xml.write('<Val></Val>\n')
            xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 9> Volts (2x1 double)
    # ----------------------------------------------------------------------

    # Pack 2x1 array of doubles.
    fmt += 'I'
    buf += struct.pack('I', 2)
    item = dict.get('VOLTS', default_volt)
    fmt += '2d'
    for i in item:
        buf += struct.pack('d', i)
    if mk_xml:
        xml.write('<Array>\n')
        xml.write('<Name>Volts</Name>\n')
        xml.write('<Dimsize>2</Dimsize>\n')
        xml.write('<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
        xml.write('</Array>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 10> Current (2x1 double)
    # ----------------------------------------------------------------------

    # Pack 2x1 array of doubles.
    fmt += 'I'
    buf += struct.pack('I', 2)
    item = dict.get('CURRENT', default_current)
    fmt += '2d'
    for i in item:
        buf += struct.pack('d', i)
    if mk_xml:
        xml.write('<Array>\n')
        xml.write('<Name>Current</Name>\n')
        xml.write('<Dimsize>2</Dimsize>\n')
        xml.write('<DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
        xml.write('</Array>\n')

    # ----------------------------------------------------------------------
    # XML Cluster closure.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('</Cluster>\n')
    return fmt, buf


def __thermal(dict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Thermal:
    # ----------------------------------------------------------------------
    default_cryostat_temp = np.zeros(8)
    default_focusbox_temp = 0

    # ----------------------------------------------------------------------
    # XML Cluster setup.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name>Thermal</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_THERMAL)
                  + '</NumElts>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 1-8> Temperature of cryostat element (double)
    # ----------------------------------------------------------------------
    # Devices are as follows:
    #   0: 70K Stage
    #   1: Hi Freq 15K Plate
    #   2: Hi Freq LNA
    #   3: Hi Freq Feedhorn
    #   4: 15K Stage
    #   5: Low Freq LNA
    #   6: 70K Radiation Shield
    #   7: Low Freq Feedhorn
    # ----------------------------------------------------------------------

    # Pack each status as unsigned int.
    item = dict.get('CRYOSTAT', default_cryostat_temp)
    for i in range(0, 8):
        fmt += 'd'
        temp_value = 0
        try:
            temp_value = item[i]
        except:
            pass
        buf += struct.pack('d', temp_value)

        # Append to XML file
        if mk_xml:
            xml.write('<DBL>\n')
            xml.write('<Name>' + THERMAL_DEF[i] + '</Name>\n')
            xml.write('<Val></Val>\n')
            xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 9> Focus Box Temperature (double)
    # ----------------------------------------------------------------------

    # Pack a double for the temperature.
    item = dict.get('FOCUSBOX', default_focusbox_temp)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>FocusBoxTemp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # XML Cluster closure.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('</Cluster>\n')
    return fmt, buf


def __receiver_lna(mydict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Receiver_LNA:
    # ----------------------------------------------------------------------
    default_lna_values = 0

    # ----------------------------------------------------------------------
    # ELEMENT 1-6> LNA Registers defined in RECEIVER_LNA_DEF (double)
    # ----------------------------------------------------------------------
    for register in RECEIVER_LNA_DEF:
        item = mydict.get(register, default_lna_values)
        buf += struct.pack('d', item)

    return fmt, buf


def __receiver(dict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Receiver:
    # ----------------------------------------------------------------------
    default_status = 0
    default_lnas = [{}, {}, {}, {}]

    # ----------------------------------------------------------------------
    # XML Cluster setup.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name>Receiver</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_RECEIVER)
                  + '</NumElts>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 1> Lo Freq Status: 0 = disable, 1 = enable (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('LOFREQSTATUS', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>LoFreqEnabled</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 2> Hi Freq Status: 0 = disable, 1 = enable (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('HIFREQSTATUS', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>HiFreqEnabled</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 3> Noise Diode Status: 0 = disable, 1 = enable (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('NOISESTATUS', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>NoiseDiodeEnabled</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 4> LNA Clusters (4x1 LNA Clusters)
    # ----------------------------------------------------------------------

    fmt += 'I'
    buf += struct.pack('I', 4)
    fmt += '12d'
    if mk_xml:
        xml.write('<Array>\n')
        xml.write('<Name>LNAs</Name>\n')
        xml.write('<Dimsize>4</Dimsize>\n')

        xml.write('<Cluster>\n')
        xml.write('<Name></Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_RECEIVER_LNA)
                  + '</NumElts>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>DrainVoltage</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>GateAVoltage</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>GateBVoltage</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>DrainCurrent</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>GateACurrent</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('<DBL>\n')
        xml.write('<Name>GateBCurrent</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

        xml.write('</Cluster>\n')
        xml.write('</Array>\n')

    item = dict.get('LNAS', default_lnas)

    # Call __receiver_lna on each LNA to generate clusters.
    for lna_cluster in item:
        append_fmt, append_buf = __receiver_lna(lna_cluster, xml, mk_xml)
        fmt += append_fmt
        buf += append_buf

    # ----------------------------------------------------------------------
    # XML Cluster closure.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('</Cluster>\n')
    return fmt, buf


def __servo_axis(dict, axis, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Servo_Axis:
    # ----------------------------------------------------------------------
    default_double = 0
    default_status = 0

    # ----------------------------------------------------------------------
    # XML Cluster setup.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name>' + str(axis) + '</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_SERVO_AXIS)
                  + '</NumElts>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 1> Positive Limit: 0 = false, 1 = true (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('POSLIMIT', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>PlusLimit</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 2> Negative Limit: 0 = false, 1 = true (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('NEGLIMIT', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>MinusLimit</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 3> Amplifier Fault: 0 = false, 1 = true (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned integer for the status.
    item = dict.get('AMPFAULT', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>AmplifierFault</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 4> Position (double)
    # ----------------------------------------------------------------------

    # Pack a double for the value.
    item = dict.get('P', default_double)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>Position</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 5> Position Error (double)
    # ----------------------------------------------------------------------

    # Pack a double for the value.
    item = dict.get('PERR', default_double)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>PositionError</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 6> Position Offset (double)
    # ----------------------------------------------------------------------

    # Pack a double for the value.
    item = dict.get('POFF', default_double)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>PositionOffset</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 7> Motor Current (double)
    # ----------------------------------------------------------------------

    # Pack a double for the value.
    item = dict.get('I', default_double)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>MotorCurrent</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # XML Cluster closure.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('</Cluster>\n')
    return fmt, buf

def __servo(dict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Servo:
    # ----------------------------------------------------------------------
    default_status = 0

    # ----------------------------------------------------------------------
    # XML Cluster setup.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('<Cluster>\n')
        xml.write('<Name>FRMServo</Name>\n')
        xml.write('<NumElts>' + str(NELEMENTS_ANTENNA_SERVO)
                  + '</NumElts>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 1> Homed: 0 = false, 1 = true (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned int for status report.
    item = dict.get('HOMED', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>Homed</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 2> SelectedRx: 0 = LF Rx, 1 = HR Rx (unsigned int)
    # ----------------------------------------------------------------------

    # Pack an unsigned int for status report.
    item = dict.get('SELECTEDRX', default_status)
    fmt += 'I'
    buf += struct.pack('I', item)
    if mk_xml:
        xml.write('<U32>\n')
        xml.write('<Name>SelectedRx</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</U32>\n')

    # ----------------------------------------------------------------------
    # ELEMENT 3-5> Axis Clusters
    # ----------------------------------------------------------------------

    # Call __servo_axis on each axis to generate clusters.
    for key, value in AXIS_DEF.items():
        item = dict.get('AXIS' + str(key), {})
        append_fmt, append_buf = __servo_axis(item, value, xml, mk_xml)
        fmt += append_fmt
        buf += append_buf

    # ----------------------------------------------------------------------
    # XML Cluster closure.
    # ----------------------------------------------------------------------
    if mk_xml:
        xml.write('</Cluster>\n')
    return fmt, buf


def __antenna(sf_dict, xml, mk_xml):
    fmt = ""
    buf = ""

    # ----------------------------------------------------------------------
    # Defaults - Antennas
    # ----------------------------------------------------------------------
    default_timestamp = 0

    sf_dict = sf_dict.get('FEM', {})
    # ----------------------------------------------------------------------
    # Dump PowerStrip
    # ----------------------------------------------------------------------
    item = sf_dict.get('POWERSTRIP', {})
    append_fmt, append_buf = __powerstrip(item, xml, mk_xml)
    fmt += append_fmt
    buf += append_buf

    # ----------------------------------------------------------------------
    # Dump Thermal
    # ----------------------------------------------------------------------
    item = sf_dict.get('THERMAL', {})
    append_fmt, append_buf = __thermal(item, xml, mk_xml)
    fmt += append_fmt
    buf += append_buf

    # ----------------------------------------------------------------------
    # Dump Receiver
    # ----------------------------------------------------------------------
    item = sf_dict.get('RECEIVER', {})
    append_fmt, append_buf = __receiver(item, xml, mk_xml)
    fmt += append_fmt
    buf += append_buf

    # ----------------------------------------------------------------------
    # Dump Servo
    # ----------------------------------------------------------------------
    item = sf_dict.get('SERVO', {})
    append_fmt, append_buf = __servo(item, xml, mk_xml)
    fmt += append_fmt
    buf += append_buf

    # ----------------------------------------------------------------------
    # Dump Timestamp. Use LabVIEW format i.e. double time in seconds since
    # 1904-01-01 00:00 UT.
    # ----------------------------------------------------------------------
    item = sf_dict.get('TIMESTAMP', default_timestamp)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>Timestamp</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    # ----------------------------------------------------------------------
    # Dump Version
    # ----------------------------------------------------------------------
    item = sf_dict.get('VERSION', VERSION)
    fmt += 'd'
    buf += struct.pack('d', item)
    if mk_xml:
        xml.write('<DBL>\n')
        xml.write('<Name>Version</Name>\n')
        xml.write('<Val></Val>\n')
        xml.write('</DBL>\n')

    return fmt, buf
