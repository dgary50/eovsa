"""
    STARBURST ACC/FEANTA Lakeshore Cryostat Worker
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""

import i_worker
import serial

# Description of Lakeshore device. Currently hard-coded.
CRYO_PORT = '/dev/ttyUSB0'
CRYO_BAUD = 9600
CRYO_BYTESIZE = 7
CRYO_STOPBITS = 1
CRYO_PARITY = 'O'
CRYO_TIMEOUT = 0.3

class CryoWorker(i_worker.IWorker):
    def __init__(self):
        super(CryoWorker, self).__init__()
        self.commands = []
        self.name = 'Cryostat-Worker'
        self.serial_connection = None

    # ---------------------------------------------------------------
    # STATEFRAME HELPERS
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: __temperature_query
        Description:
            Queries the temperatures monitored by the Lakeshore and
            returns them in a single, parsed array.
        Returns:
            returnVal: array of floats representing the temperatures in
                degrees Kelvin.
    """
    # endregion
    def __temperature_query(self):
        query_cmd = 'krdg? 0\x0d\x0a'
        self.serial_connection = serial.Serial(
            port=CRYO_PORT, baudrate=CRYO_BAUD, bytesize=CRYO_BYTESIZE,
            parity=CRYO_PARITY, stopbits=CRYO_STOPBITS, timeout=CRYO_TIMEOUT)
        self.serial_connection.write(query_cmd)
        returnString = self.serial_connection.readline()
        returnString = returnString.split(',')
        returnVal = []
        for number in returnString:
            try:
                returnVal.append(float(number))
            except ValueError:
                returnVal.append(0)
        self.serial_connection.close()
        self.serial_connection = None
        return returnVal

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
        pass

    # region Method Description
    """
    Method: stateframe_query
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def stateframe_query(self):
        temperatures = self.__temperature_query()
        return {'CRYOSTAT': temperatures}