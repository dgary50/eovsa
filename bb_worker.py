"""
    STARBURST ACC/FEANTA BeagleBone Worker
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""

import i_worker
import numpy as np
import socket

# Description of the BeagleBone device. Currently hard-coded.
BB_HOSTNAME = 'lna14.solar.pvt'
BB_PORT = 50002
BB_TIMEOUT = 0.3

# Scale Factors
DRAIN_FACTOR = 0.300
GATE_FACTOR = 0.388
CURRENT_FACTOR = -2

# Query Dictionary
QUERY_DICT = {0: 'DRAINVOLTAGE',
              2: 'GATEAVOLTAGE',
              4: 'GATEBVOLTAGE',
              1: 'DRAINCURRENT',
              3: 'GATEACURRENT',
              5: 'GATEBCURRENT'}

# Amp-number Mapping
AMP_MAP = {'hh': 0,
           'lh': 1,
           'lv': 2,
           'hv': 3}

class BBWorker(i_worker.IWorker):
    def __init__(self):
        super(BBWorker, self).__init__()
        self.commands = ['LNA-GATE1',
                         'LNA-GATE2',
                         'LNA-DRAIN',
                         'LNA-ENABLE']
        self.name = 'BB-Worker'
        self.bb_socket = None
        self.bb_ip = socket.gethostbyname(BB_HOSTNAME)
        self.dt = np.dtype('float32').newbyteorder('>')

    # ---------------------------------------------------------------
    # COMMAND ROUTINES
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: __lna_gate1
        Description:
            Routine to build command to change voltage on the LNA
            Gate A.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['LNA-GATE1', amp_number, voltage]
        Returns:
            command: command designated to complete task.
    """
    # endregion
    def __lna_gate1(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 3:
            self.logger('Invalid call to LNA-GATE1.')
            return None
        amp_num = acc_command[1]
        voltage = None
        try:
            amp_num = AMP_MAP.get(amp_num.lower(), -1)
            voltage = float(acc_command[2])
            voltage /= GATE_FACTOR
        except ValueError:
            self.logger('Invalid call to LNA-GATE1.')
            return None
        if amp_num == -1:
            self.logger('Invalid call to LNA-GATE1.')
            return None

        # Given that the parameters are all reasonable, we return the
        # command string to be processed later.
        command = ['set amp ' + str(amp_num) + ' gatea ' + str(voltage),
                   'latch']
        return command

    # region Method Description
    """
    Method: __lna_gate2
        Description:
            Routine to build command to change voltage on the LNA
            Gate B.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['LNA-GATE2', amp_number, voltage]
        Returns:
            command: command designated to complete task.
    """
    # endregion
    def __lna_gate2(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 3:
            self.logger('Invalid call to LNA-GATE2.')
            return None
        amp_num = acc_command[1]
        voltage = None
        try:
            amp_num = AMP_MAP.get(amp_num.lower(), -1)
            voltage = float(acc_command[2])
            voltage /= GATE_FACTOR
        except ValueError:
            self.logger('Invalid call to LNA-GATE2.')
            return None
        if amp_num == -1:
            self.logger('Invalid call to LNA-GATE2.')
            return None

        # Given that the parameters are all reasonable, we return the
        # command string to be processed later.
        command = ['set amp ' + str(amp_num) + ' gateb ' + str(voltage),
                   'latch']
        return command

    # region Method Description
    """
    Method: __lna_drain
        Description:
            Routine to build command to change voltage on the LNA drain.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['LNA-DRAIN', amp_number, voltage]
        Returns:
            command: command designated to complete task.
    """
    # endregion
    def __lna_drain(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 3:
            self.logger('Invalid call to LNA-DRAIN.')
            return None
        amp_num = acc_command[1]
        voltage = None
        try:
            amp_num = AMP_MAP.get(amp_num.lower(), -1)
            voltage = float(acc_command[2])
            voltage /= DRAIN_FACTOR
        except ValueError:
            self.logger('Invalid call to LNA-DRAIN.')
            return None
        if amp_num == -1:
            self.logger('Invalid call to LNA-DRAIN.')
            return None

        # Given that the parameters are all reasonable, we return the
        # command string to be processed later.
        command = ['set amp ' + str(amp_num) + ' drain ' + str(voltage),
                   'latch']
        return command

    # region Method Description
    """
    Method: __lna_enable
        Description:
            Routine to set power to LNA.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['LNA-ENABLE', amp_number, state]
        Returns:
            command: command designated to complete task.
    """
    # endregion
    def __lna_enable(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 3:
            self.logger('Invalid call to LNA-ENABLE.')
            return None
        amp_num = acc_command[1]
        state = acc_command[2]
        try:
            amp_num = AMP_MAP.get(amp_num.lower(), -1)
        except ValueError:
            self.logger('Invalid call to LNA-ENABLE.')
            return None
        if amp_num == -1:
            self.logger('Invalid call to LNA-ENABLE.')
            return None
        if state.lower() == 'on':
            state = 1
        elif state.lower() == 'off':
            state = 0
        else:
            self.logger('Invalid call to LNA-ENABLE.')
            return None

        # Given that the parameters are all reasonable, we return the
        # command string to be processed later.
        command = ['set power ' + str(amp_num) + ' ' + str(state)]
        return command

    # ---------------------------------------------------------------
    # FUNCTION MAP
    # ---------------------------------------------------------------
    function_map = {'LNA-GATE1': __lna_gate1,
                    'LNA-GATE2': __lna_gate2,
                    'LNA-DRAIN': __lna_drain,
                    'LNA-ENABLE': __lna_enable}

    # ---------------------------------------------------------------
    # STATEFRAME HELPERS
    # ---------------------------------------------------------------
    def __lna_query(self):
        query_cmd = 'read\r\n'
        query_socket = socket.socket(socket.AF_INET,
                                     socket.SOCK_STREAM)
        query_socket.settimeout(BB_TIMEOUT)
        query_socket.connect((self.bb_ip, BB_PORT))
        query_socket.sendall(query_cmd)
        read_buf = query_socket.recv(96)
        query_socket.close()
        data = np.fromstring(read_buf, self.dt)

        amp0 = {}
        amp1 = {}
        amp2 = {}
        amp3 = {}

        for i in range(0, 6):
            scale = 1
            if i % 2 == 1:
                scale = CURRENT_FACTOR
            amp0[QUERY_DICT[i]] = data[i * 4] / scale
            amp1[QUERY_DICT[i]] = data[i * 4 + 1] / scale
            amp2[QUERY_DICT[i]] = data[i * 4 + 2] / scale
            amp3[QUERY_DICT[i]] = data[i * 4 + 3] / scale

        return [amp0, amp1, amp2, amp3]

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

        # Use the routine calls to generate url commands.
        command_strings = self.function_map[acc_command[0]](self, acc_command)
        if command_strings is not None:
            for command_string in command_strings:
                self.bb_socket = socket.socket(socket.AF_INET,
                                               socket.SOCK_STREAM)
                self.bb_socket.settimeout(BB_TIMEOUT)
                self.bb_socket.connect((self.bb_ip, BB_PORT))
                self.logger('The following command was issued: ' +
                            command_string)

                command_string += '\r\n'
                self.bb_socket.sendall(command_string)
                self.bb_socket.close()
                self.bb_socket = None

    # region Method Description
    """
    Method: stateframe_query
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def stateframe_query(self):
        lnas = self.__lna_query()
        return {'LNAS': lnas}