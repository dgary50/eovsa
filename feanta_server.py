"""
    STARBURST ACC/FEANTA Middle Server
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""

import datetime
import socket
import sys
import time
import threading
import gen_fem_sf
import traceback

# Logging information.
TIMESTAMP_FMT = '%Y-%m-%d %H:%M:%S'
LOG_FILE = 'bridge_server.log'

# Define all constants:
# Currently hard-coded, will eventually be read from acc.ini
HOST = ''
HOST_PORT = 5676
ACC_HOSTNAME = 'acc.solar.pvt'
ACC_PORT = 5675
VERSION = 1.2  # Version date: 10/6/2015


# region Class Description
"""
Class: ServerDaemon
    Description:
        Implementation of a daemon server that is intended to run on
        the feanta box in order to process commands from the ACC and
        direct and execute the commands to the sub-units connected to
        the feanta computer.
    Arguments:
        pidfile: string designating the .pid file to save the pid for
            this daemon process to allow for the process to be stopped
            by the stop function or to be stopped manually in Linux.
"""
# endregion
class ServerDaemon():
    def __init__(self, pidfile):
        self.pidfile_path = pidfile
        self.pidfile_timeout = 5
        self.stdin_path = '/dev/null'
        self.stdout_path = '/dev/null'
        self.stderr_path = '/dev/null'
        self.workers = {}
        self.function_map = {}
        self.log_file = LOG_FILE
        self.acc_ip = socket.gethostbyname(ACC_HOSTNAME)

    # ---------------------------------------------------------------
    # BASIC ROUTINES:
    # ---------------------------------------------------------------
    def __get_timestamp(self):
        current_time = time.time()
        timestamp = datetime.datetime.fromtimestamp(current_time)
        timestamp = timestamp.strftime(TIMESTAMP_FMT)
        return timestamp

    def __log(self, message):
        log_message = self.__get_timestamp() + ': ' + str(message) + '\n'
        f = open(self.log_file, "a")
        f.write(log_message)
        f.close()
        print log_message

    # ---------------------------------------------------------------
    # CORE ROUTINES
    # ---------------------------------------------------------------

    # region Method Description
    """
    Method: link_worker
        Description:
            This method is used to link a worker extending the i_worker
            class to this server so that commands associated with the
            i_worker can be executed properly through this server.
        Arguments:
            worker: the target worker to be linked to this server
    """
    # endregion
    def link_worker(self, worker):
        self.workers[worker.name] = worker
        for command in worker.get_command_list():
            self.function_map[command] = worker
        worker.set_logger(self.__log)

    # region Method Description
    """
    Method: list_workers
        Description:
            Lists each worker linked to this ServerDaemon.
    """
    # endregion
    def list_workers(self):
        workers_list = ''
        for worker in self.workers.keys():
            workers_list += worker + '\n'
        return workers_list

    # region Method Description
    """
    Method: set_log_file
        Description:
            Sets the destination path for the log file. Defaulted to
            LOG_FILE.
    """
    # endregion
    def set_log_file(self, log_file_destination):
        self.log_file = log_file_destination

    # region Method Description
    """
    Method: list_commands
        Description:
            Lists every command that this server can respond to.
    """
    # endregion
    def list_commands(self):
        return self.function_map.keys()

    def make_stateframe_dict(self):
        fem_dict = {}

        # Handle powerstrip cluster.
        worker = self.workers.get('PDU-Worker', None)
        if worker is not None:
            try:
                fem_dict['POWERSTRIP'] = worker.stateframe_query()
            except Exception, e:
                fem_dict['POWERSTRIP'] = {}
                self.__log(traceback.format_exc())
        else:
            fem_dict['POWERSTRIP'] = {}

        # Handle thermal cluster.
        worker = self.workers.get('Cryostat-Worker', None)
        working_dict = {}
        if worker is not None:
            try:
                working_dict = worker.stateframe_query()
            except Exception, s:
                self.__log(traceback.format_exc())

        worker = self.workers.get('Temp-Worker', None)
        if worker is not None:
            try:
                working_dict['FOCUSBOX'] = worker.stateframe_query()
            except Exception, e:
                working_dict['FOCUSBOX'] = 0
                self.__log(traceback.format_exc())
        else:
            working_dict['FOCUSBOX'] = 0

        fem_dict['THERMAL'] = working_dict

        # Handle receiver cluster.
        worker = self.workers.get('BB-Worker', None)
        working_dict = {}
        if worker is not None:
            try:
                working_dict = worker.stateframe_query()
            except Exception, e:
                pass

        working_dict['LOFREQSTATUS'] = 0
        working_dict['HIFREQSTATUS'] = 0
        working_dict['NOISESTATUS'] = 0
        fem_dict['RECEIVER'] = working_dict

        # Handle servo cluster.
        worker = self.workers.get('GeoBrick-Worker', None)
        if worker is not None:
            try:
                fem_dict['SERVO'] = worker.stateframe_query()
            except Exception, e:
                fem_dict['SERVO'] = {}
                self.__log(traceback.format_exc())
        else:
            fem_dict['SERVO'] = {}

        # Handle version.
        fem_dict['VERSION'] = VERSION

        # Handle timestamp
        fem_dict['TIMESTAMP'] = time.time() + 2082844800

        return {'FEM': fem_dict}

    def send_stateframe_dict(self):
        try:
            fem_dict = self.make_stateframe_dict()
            fmt, buf, xml = gen_fem_sf.gen_fem_sf(fem_dict)
            packet_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            packet_socket.settimeout(0.3)
            packet_socket.connect((self.acc_ip, ACC_PORT))
            packet_socket.sendall(buf)
            packet_socket.close()
            # persec = open('/tmp/persec.txt', 'a')
            # persec.write(buf + '\n')
            # persec.close()
        finally:
            threading.Timer(0.3, self.send_stateframe_dict).start()


    # region Method Description
    """
    Method: run
        Description:
            Daemon routine for the ServerDaemon between the ACC and the Brick.
    """
    # endregion
    def run(self):
        # Setup listener to this box at HOST_PORT.
        acc_listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.__log('Attempt to set up listener...')
        try:
            acc_listener.bind((HOST, HOST_PORT))
        except socket.error, msg:
            self.__log('Unable to listen at port ' + str(HOST_PORT) +
                       '. Error Code: ' + str(msg[0]) + '. Message: ' +
                       str(msg[1]))
            sys.exit()
        acc_listener.listen(1)
        self.__log('Successfully setup listener')

        polling_thread = threading.Thread(target=self.send_stateframe_dict)
        polling_thread.start()

        while True:
            # Wait for a connection from ACC.
            connection, address = acc_listener.accept()
            self.__log('Connection from ' + address[0] +
                       ':' + str(address[1]))

            # Read packet sent from ACC, currently capped at 1024 bytes.
            acc_command = connection.recv(1024)
            self.__log('Command issued from connection: ' + acc_command)
            acc_command = acc_command.split()

            # Echo command issued back to ACC.
            connection.sendall(acc_command[0])

            # Verify that the given command exists and execute it with the
            # correct worker if it does.
            try:
                worker = self.function_map[acc_command[0]]
                worker.execute(acc_command)
            except KeyError:
                self.__log('Unrecognized command received: ' +
                           acc_command[0] + '.')
