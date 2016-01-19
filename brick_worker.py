"""
    STARBURST ACC/FEANTA GeoBrick Worker
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""
# Change-log:
#   2015-12-30  DG
#     Changed __frm_rx_sel() to use RX values 'LO' and 'HI' 
#     instead of 1, 2.  Also eliminated unhandled error when
#     brick returns bad information.  Also changed RXSEL to
#     SELECTEDRX to reflect change in gen_fem_sf.py
#   2016-01-01  DG
#     Trap unhandled timeout exception in reading from BRICK
#     (in case BRICK outlet is powered off). Added try..except clause 
#     to just return zeros in __brickmonitor_query(), for this case.
#
import i_worker
import socket
import struct

# Description of the GeoBrick device. Currently hard-coded.
BRICK_HOSTNAME = 'geobrickanta.solar.pvt'
BRICK_PORT = 1025
BRICK_TIMEOUT = 0.5

# Program spaces that can be used in the GeoBrick.
COMMAND_REGIS = 'P1000='
ARG1_REGIS = ' P1001='
ARG2_REGIS = ' P1002='

# Brick command dictionary.
COMMAND_DICT = {'Home': 1,
                'SelRx': 2,
                'SetAngle': 3,
                'SetZOffset': 4,
                'SetXOffset': 5,
                'Kill': 6,
                'Enable': 7,
                'SetX': 8,
                'SetZ': 9}

# Dictionaries for ethernet packets to the Brick.
RQ_TYPE = {'upload': '\xc0',
           'download': '\x40'}
RQ = {'sendline': '\xb0',
      'getline': '\xb1',
      'flush': '\xb3',
      'getmem': '\xb4',
      'setmem': '\xb5',
      'setbit': '\xba',
      'setbits': '\xbb',
      'port': '\xbe',
      'getresponse': '\xbf',
      'readready': '\xc2',
      'response': '\xc4',
      'getbuffer': '\xc5',
      'writebuffer': '\xc6',
      'writeerror': '\xc7',
      'fwdownload': '\xcb',
      'ipaddress': '\xe0'}
COORDINATE = {1: 'Z',
              3: 'A',
              4: 'X'}
RX = {'LO':1,
      'HI':2}
AXIS_SCALING = {1: 42.5636 * 96 * 32,
                3: 23181.5208 * 96 * 32,
                4: 3973.477 * 96 * 32}
MPADDRESSSTART = 900

class BrickWorker(i_worker.IWorker):
    def __init__(self):
        super(BrickWorker, self).__init__()
        self.commands = ['FRM-HOME',
                         'FRM-KILL',
                         'FRM-RX-SEL',
                         'FRM-SET-PA',
                         'FRM-X-OFFSET',
                         'FRM-Z-OFFSET',
                         'FRM-ABS-X',
                         'FRM-ABS-Z',
                         'FRM-ENABLE']
        self.brick_socket = None
        self.brick_ip = socket.gethostbyname(BRICK_HOSTNAME)
        self.name = 'GeoBrick-Worker'

    # ---------------------------------------------------------------
    # COMMAND PACKAGING ROUTINES SPECIFIC TO GEOBRICK
    # ---------------------------------------------------------------

    #region Method Description
    """
    Method: __make_brick_command
        Description:
            Takes a command to the Brick and packages it into an
            ethernet packet recognized by the Brick system.
        Arguments:
            rq_type: type of request, either 'upload' or 'download'.
            rq: nature of request, lookup dictionary defined in RQ.
            val: value associated with the request.
            index: index associated with the request.
            command_packets: list of strings to be packed into TCP packets.
    """
    #endregion
    def __make_brick_command(self, rq_type, rq, val, index, command_packets):
        packets = []
        for packet in command_packets:
            buf = RQ_TYPE[rq_type] + RQ[rq]
            buf += struct.pack('H', val)
            buf += struct.pack('H', index)
            buf += struct.pack('H', socket.htons(len(packet) + 1))
            buf += struct.pack(str(len(packet)) + 's', packet)
            buf += struct.pack("B", 0)
            packets.append(buf)
        return packets

    # ---------------------------------------------------------------
    # COMMAND ROUTINES
    # ---------------------------------------------------------------

    #region Method Description
    """
    Method: __frm_home
        Description:
            Runs homing procedure local to the GeoBrick. Do NOT use this
            method on its own. This method is error checked before execution.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['FRM-HOME']
        Returns:
            [0]: A list of packets as strings before compression.
            [1]: A list of TCP/Ethernet packets ready to be sent to the Brick.
    """
    #endregion
    def __frm_home(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 1:
            self.logger('Invalid call to FRM-HOME.')
            return None

        command_packets = []

        command = COMMAND_REGIS + str(COMMAND_DICT['Home'])
        command_packets.append(command)

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                         0, 0, command_packets)

    #region Method Description
    """
    Method: __frm_rx_sel
        Description:
            Routine to select one of two receivers on the antenna
        Arguments:
            acc_command: list of strings sent from the ACC. List format:
                ['FRM-RX-SEL', rx] where rx is 1 for low-nu and 2 for high-nu.
        Returns:
            [0]: A list of packets as strings before compression.
            [1]: A list of TCP/Ethernet packets ready to be sent to the Brick.
    """
    #endregion
    def __frm_rx_sel(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-RX-SEL.')
            return None
        rxstr = None
        try:
            rxstr = acc_command[1].upper()
            if rxstr not in ['LO', 'HI']:
                raise ValueError('Invalid RX selection: '+rxstr)
        except ValueError:
            self.logger('Invalid call to FRM-RX-SEL.')
            return None

        # Build command based on parameters.
        command = COMMAND_REGIS + str(COMMAND_DICT['SelRx']) + \
                  ARG1_REGIS + str(RX[rxstr])
        command_packets = [command]

        return command_packets, self.__make_brick_command('download',
                                                          'getresponse',
                                                          0, 0,
                                                          command_packets)

    #region Method Description
    """
    Method: __frm_set_pa
        Description:
            Routine to move motor 3 to a given angle. This routine should only
            be called after a FRM_HOME command has been issued. Do NOT use
            this method on its own.
        Arguments:
            acc_command: list of strings sent from the ACC. List format:
                ['FRM-SET-PA', angle] where angle is the absolute angle to be
                set.
        Returns:
            [0]: A list of packets as strings before compression.
            [1]: A list of TCP/Ethernet packets ready to be sent to the Brick.
    """
    #endregion
    def __frm_set_pa(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-SET-PA.')
            return None
        angle = None
        try:
            angle = int(acc_command[1])
            if angle > 90 or angle < -90:
                raise ValueError('Invalid position angle selection.')
        except ValueError:
            self.logger('Invalid call to FRM-SET-PA.')
            return None

        # Build command based on parameters.
        command = COMMAND_REGIS + str(COMMAND_DICT['SetAngle']) + \
                  ARG1_REGIS + str(angle)
        command_packets = [command]

        return command_packets, self.__make_brick_command('download',
                                                          'getresponse',
                                                          0, 0,
                                                          command_packets)

    def __frm_x_offset(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-X-OFFSET.')
            return None
        offset = None
        try:
            offset = float(acc_command[1])
        except ValueError:
            self.logger('Invalid call to FRM-X-OFFSET.')
            return None

        command = COMMAND_REGIS + str(COMMAND_DICT['SetXOffset']) + \
                      ARG1_REGIS + str(offset)

        # Build command based on parameters. (This assumes that the
        # position given is in physical units.)
        command_packets = [command]

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                         0, 0, command_packets)

    def __frm_z_offset(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-Z-OFFSET.')
            return None
        offset = None
        try:
            offset = float(acc_command[1])
        except ValueError:
            self.logger('Invalid call to FRM-Z-OFFSET.')
            return None

        command = COMMAND_REGIS + str(COMMAND_DICT['SetZOffset']) + \
                      ARG1_REGIS + str(offset)

        # Build command based on parameters. (This assumes that the
        # position given is in physical units.)
        command_packets = [command]

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                         0, 0, command_packets)

    #region Method Description
    """
    Method: __frm_abs_x
        Description:
            Routine to move x-axis to specified location
        Arguments:
            acc_command: list of strings sent from the ACC. List format:
                ['FRM-ABS-X', destination] where destination is the destination
                in physical units (mm).
        Returns:
            [0]: A list of packets as strings before compression.
            [1]: A list of TCP/Ethernet packets ready to be sent to the Brick.
    """
    #endregion
    def __frm_abs_x(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-ABS-X.')
            return None
        position = None
        try:
            position = float(acc_command[1])
        except ValueError:
            self.logger('Invalid call to FRM-ABS-X.')
            return None

        command = COMMAND_REGIS + str(COMMAND_DICT['SetX']) + \
                      ARG1_REGIS + str(position)

        # Build command based on parameters. (This assumes that the
        # position given is in physical units.)
        command_packets = [command]

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                         0, 0, command_packets)

    #region Method Description
    """
    Method: __frm_abs_z
        Description:
            Routine to move z-axis to specified location
        Arguments:
            acc_command: list of strings sent from the ACC. List format:
                ['FRM-ABS-Z', destination] where destination is the destination
                in physical units (mm).
        Returns:
            [0]: A list of packets as strings before compression.
            [1]: A list of TCP/Ethernet packets ready to be sent to the Brick.
    """
    #endregion
    def __frm_abs_z(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 2:
            self.logger('Invalid call to FRM-ABS-Z.')
            return None
        position = None
        try:
            position = float(acc_command[1])
        except ValueError:
            self.logger('Invalid call to FRM-ABS-Z.')
            return None

        command = COMMAND_REGIS + str(COMMAND_DICT['SetZ']) + \
                      ARG1_REGIS + str(position)

        # Build command based on parameters. (This assumes that the
        # position given is in physical units.)
        command_packets = [command]

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                         0, 0, command_packets)

    def __frm_kill(self, acc_command):
        # Error check that the command given was formatted correctly.
        if len(acc_command) != 1:
            self.logger('Invalid call to FRM-KILL')
            return None

        command_packets = []

        command = COMMAND_REGIS + str(COMMAND_DICT['Kill'])
        command_packets.append(command)

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                        0, 0, command_packets)

    def __frm_enable(self, acc_command):
        # Error check that the command given was formatted correctly.
        if len(acc_command) != 1:
            self.logger('Invalid call to FRM-ENABLE')
            return None

        command_packets = []

        command = COMMAND_REGIS + str(COMMAND_DICT['Enable'])
        command_packets.append(command)

        return command_packets, \
               self.__make_brick_command('download', 'getresponse',
                                        0, 0, command_packets)

    # ---------------------------------------------------------------
    # FUNCTION MAP
    # ---------------------------------------------------------------
    function_map = {'FRM-HOME': __frm_home,
                    'FRM-KILL': __frm_kill,
                    'FRM-RX-SEL': __frm_rx_sel,
                    'FRM-SET-PA': __frm_set_pa,
                    'FRM-X-OFFSET': __frm_x_offset,
                    'FRM-Z-OFFSET': __frm_z_offset,
                    'FRM-ABS-X': __frm_abs_x,
                    'FRM-ABS-Z': __frm_abs_z,
                    'FRM-ENABLE': __frm_enable}

    # ---------------------------------------------------------------
    # STATEFRAME HELPERS
    # ---------------------------------------------------------------
    def __brickmonitor_query(self):
        command = 'LIST GATHER'
        query_socket = socket.socket(socket.AF_INET,
                                     socket.SOCK_STREAM)
        query_socket.settimeout(BRICK_TIMEOUT)
        try:
            query_socket.connect((self.brick_ip, BRICK_PORT))
        except socket.timeout:
            # Connect times out when, for example, the BRICK outlet
            # is powered off.
            self.logger(
                'Unable to connect to BRICK--returning zeros')
            return [0.0]
        cmd_string = [command]
        cmd = self.__make_brick_command('download', 'getresponse',
                                        0, 0, cmd_string)
        query_socket.sendall(cmd[0])
        response = query_socket.recv(1024)
        query_socket.close()
        response = response.replace('\r', ' ')
        response = response.split(' ')
        parsed_response = []
        for monitor_point in response:
            parsed_response.append(self.__str2float(monitor_point))

        return parsed_response, response

    # Expose direct information from Brick, for debugging
    def brickmonitor_query(self):
        parsed_response, response = self.__brickmonitor_query()
        return parsed_response, response
        
    def __str2float(self, str_val):
        num = 0
        try:
            num = int(str_val, 16)
        except Exception:
            num = 0
        return (num >> 12) * 2**((num & 0xFFF) - 2082)

    # ---------------------------------------------------------------
    # INTERFACE IMPLEMENTATIONS
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: get_command_list
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def get_command_list(self):
        return self.commands

    # region Method Description
    """
    Method: execute
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def execute(self, acc_command):
        # Use the routine functions to get the commands to push.
        packets = self.function_map[acc_command[0]](
                    self, acc_command)

        if packets is not None:
            self.logger('Issued the following commands to brick:')
            for packet in packets[0]:
                self.logger(repr(packet))

            # Try pushing message across TCP.
            # Wait for reply of at most 1024 bytes.
            try:
                for packet in packets[1]:
                    reply = None
                    self.brick_socket = socket.socket(socket.AF_INET,
                                              socket.SOCK_STREAM)
                    self.brick_socket.connect((self.brick_ip, BRICK_PORT))
                    self.brick_socket.sendall(packet)
                    self.brick_socket.settimeout(BRICK_TIMEOUT)
                    reply = self.brick_socket.recv(1024)

                    self.logger('Reply from brick: ' + reply)
                    self.brick_socket.close()
                    self.brick_socket = None
            except socket.gaierror:
                self.logger('Brick hostname could not be resolved.')
            except socket.error:
                self.logger('Unable to send packet to brick.')

    # region Method Description
    """
    Method: stateframe_query
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def stateframe_query(self):
        stateframe_data = {'AXIS1': {},
                           'AXIS3': {},
                           'AXIS4': {}}
        fetched_data, junk = self.__brickmonitor_query()
        # Check if fetched_data is truncated, and zero-fill if so
        if len(fetched_data) < 24:
            fetched_data += [0.0]*(24-len(fetched_data))
        stateframe_data['HOMED'] = \
            int(fetched_data[1])
        stateframe_data['SELECTEDRX'] = \
            int(fetched_data[2])

        stateframe_data['AXIS1']['P'] = \
            float(fetched_data[3])
        stateframe_data['AXIS1']['PERR'] = \
            float(fetched_data[4])
        stateframe_data['AXIS1']['POFF'] = \
            float(fetched_data[5])
        stateframe_data['AXIS1']['I'] = \
            float(fetched_data[6])
        stateframe_data['AXIS1']['POSLIMIT'] = \
            int(fetched_data[7])
        stateframe_data['AXIS1']['NEGLIMIT'] = \
            int(fetched_data[8])
        stateframe_data['AXIS1']['AMPFAULT'] = \
            int(fetched_data[9])

        stateframe_data['AXIS3']['P'] = \
            float(fetched_data[10])
        stateframe_data['AXIS3']['PERR'] = \
            float(fetched_data[11])
        stateframe_data['AXIS3']['POFF'] = \
            float(fetched_data[12])
        stateframe_data['AXIS3']['I'] = \
            float(fetched_data[13])
        stateframe_data['AXIS3']['POSLIMIT'] = \
            int(fetched_data[14])
        stateframe_data['AXIS3']['NEGLIMIT'] = \
            int(fetched_data[15])
        stateframe_data['AXIS3']['AMPFAULT'] = \
            int(fetched_data[16])

        stateframe_data['AXIS4']['P'] = \
            float(fetched_data[17])
        stateframe_data['AXIS4']['PERR'] = \
            float(fetched_data[18])
        stateframe_data['AXIS4']['POFF'] = \
            float(fetched_data[19])
        stateframe_data['AXIS4']['I'] = \
            float(fetched_data[20])
        stateframe_data['AXIS4']['POSLIMIT'] = \
            int(fetched_data[21])
        stateframe_data['AXIS4']['NEGLIMIT'] = \
            int(fetched_data[22])
        stateframe_data['AXIS4']['AMPFAULT'] = \
            int(fetched_data[23])

        return stateframe_data
