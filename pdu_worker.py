"""
    STARBURST ACC/FEANTA PDU Worker
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""
# Change-log:
#   2016-01-01  DG
#     During the rare sending of a command to the PDU, the stateframe
#     query was failing.  Added try..except clause to just return zeros
#     in __statusandpower_query(), for this case.
#
from bs4 import BeautifulSoup as Soup
import mechanize
import urllib
import i_worker
import threading

# Description of the PDU device. Currently hard-coded.
PDU_HOSTNAME = 'http://pduanta.solar.pvt'
PDU_USERNAME = 'admin'
PDU_PASSWORD = 'pwr4me'

ON_OFF_MAP = {0: 'OFF',
              1: 'ON',
              'ON': 1,
              'OFF': 0}

LOGIN_DATA = {'Username': PDU_USERNAME,
              'Password': PDU_PASSWORD}

PDU_TIMEOUT = 0.3


class PDUWorker(i_worker.IWorker):
    def __init__(self):
        super(PDUWorker, self).__init__()
        self.commands = ['OUTLET',
                         'ND-ON',
                         'ND-OFF']
        self.browser = None
        self.name = 'PDU-Worker'
        self.lock = threading.Lock()

    # ---------------------------------------------------------------
    # LOGIN ROUTINES SPECIFIC TO PDU
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: __login
        Description:
            Method used to login to the PDU using the parameters
            defined at the top of this file.
    """
    # endregion
    def __login(self):
        try:
            self.browser.open(PDU_HOSTNAME + '/index.htm', timeout=PDU_TIMEOUT)
            if self.browser.geturl() == PDU_HOSTNAME + '/index.htm':
                return True
            else:
                raise Exception()
        except:
            self.browser = mechanize.Browser()
            self.browser.set_handle_robots(False)
            self.browser.set_handle_refresh(False)

            encoded_data = urllib.urlencode(LOGIN_DATA)
            self.browser.open(PDU_HOSTNAME + '/login.tgi',
                              encoded_data, timeout=PDU_TIMEOUT)
            self.browser.open(PDU_HOSTNAME + '/index.htm', timeout=PDU_TIMEOUT)
            if self.browser.geturl() == PDU_HOSTNAME + '/index.htm':
                self.logger('Successfully logged into PDU.')
                return True
        return False

    # region Method Description
    """
    Method: __logout
        Description:
            Method used to logout of a connection to the PDU after execution
            of procedures.
    """
    # endregion
    def __logout(self):
        self.browser.open(PDU_HOSTNAME + '/logout')
        self.browser.close()
        self.browser = None
        self.logger('Successfully logged out.')

    # ---------------------------------------------------------------
    # COMMAND ROUTINES
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: __outlet
        Description:
            Routine to build url that will switch a designated outlet
            on the PDU on/off.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['OUTLET', outlet_number, 'on' or 'off']
        Returns:
            command: url designated to complete task.
    """
    # endregion
    def __outlet(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 3:
            self.logger('Invalid call to OUTLET.')
            return None
        outlet_num = None
        on_off = acc_command[2]
        try:
            outlet_num = int(acc_command[1])
            if (outlet_num < 1) or (outlet_num > 8) or \
                    (on_off.upper() != 'ON' and on_off.upper() != 'OFF'):
                self.logger('Invalid call to OUTLET.')
                return None
        except ValueError:
            self.logger('Invalid call to OUTLET.')
            return None

        # Given that the parameters are all correct, we return the
        # link string to be processed later.
        command = PDU_HOSTNAME + '/outlet?' + str(outlet_num) + '=' + \
                  on_off.upper()
        return command

    # region Method Description
    """
    Method: __nd_on
        Description:
            Routine to build url that specifically switches the noise diode
            (outlet number 8) on.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['ND-ON']
        Returns:
            command: url designated to complete task.
    """
    # endregion
    def __nd_on(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 1:
            self.logger('Invalid call to ND-ON.')
            return None
        # Given that the parameters are all correct, we return the
        # link string to be processed later.
        command = PDU_HOSTNAME + '/outlet?' + str(8) + '=ON'
        return command

    # region Method Description
    """
    Method: __nd_off
        Description:
            Routine to build url that specifically switches the noise diode
            (outlet number 8) off.
        Arguments:
            acc_command: list of the strings sent from the ACC. List format:
                ['ND-OFF']
        Returns:
            command: url designated to complete task.
    """
    # endregion
    def __nd_off(self, acc_command):
        # Error check that the command given is formatted correctly.
        if len(acc_command) != 1:
            self.logger('Invalid call to ND-OFF.')
            return None
        # Given that the parameters are all correct, we return the
        # link string to be processed later.
        command = PDU_HOSTNAME + '/outlet?' + str(8) + '=OFF'
        return command

    # ---------------------------------------------------------------
    # FUNCTION MAP
    # ---------------------------------------------------------------
    function_map = {'OUTLET': __outlet,
                    'ND-ON': __nd_on,
                    'ND-OFF': __nd_off}

    # ---------------------------------------------------------------
    # STATEFRAME HELPERS
    # ---------------------------------------------------------------
    def __statusandpower_query(self):
        statuses = []
        volts = []
        current = []
        is_logged_in = self.__login()
        if not is_logged_in:
            self.logger('Unable to login to PDU.')

        else:
            # Get statuses of the 8 devices.
            html_response = self.browser.response()
            try:
                xml_data = html_response.read()
            except mechanize.HTTPError:
                # Read fails when trying to read while sending a
                # command, so just return zeros in that case.
                self.logger(
                    'PDU statusandpower query failed--returning zeros')
                return [0,0,0,0,0,0,0,0], [0.0,0.0], [0.0,0.0]
            xml_soup = Soup(xml_data, 'html.parser')
            read_statuses = xml_soup('table')[5]('font')
            for status in read_statuses:
                statuses.append(ON_OFF_MAP[status.text])

            # Get voltage and current of the power strip
            read_power = xml_soup('table')[5]('th', {'colspan': '3'})
            for power_reading in read_power:
                readings = power_reading.text.split()
                v = 0
                i = 0
                try:
                    v = float(readings[0].replace('V', ''))
                    i = float(readings[1].replace('A', ''))
                except ValueError:
                    pass
                volts.append(v)
                current.append(i)

        return statuses, volts, current

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
        is_logged_in = self.__login()
        if is_logged_in:

            # Use the routine calls to generate url commands.
            command_string = self.function_map[acc_command[0]]\
                                (self, acc_command)
            if command_string is not None:
                self.browser.open(command_string)
                self.logger(
                    'The following link was followed for the PDU: ' +
                    command_string)
        else:
            self.logger('Unable to login to PDU.')

    # region Method Description
    """
    Method: stateframe_query
        Description:
            Refer to abstract class IWorker located in i_worker.py
            for full description.
    """
    # endregion
    def stateframe_query(self):
        statuses, volt, current = self.__statusandpower_query()
        return {'STATUS': statuses,
                'VOLTS': volt,
                'CURRENT': current}
