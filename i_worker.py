"""
    STARBURST ACC/FEANTA Worker Interface
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""

# region Class Description
"""
Class: IWorker
    Description:
        IWorker is an abstract class that defines the necessary methods for
        a worker that can be used by the ServerDaemon class in feanta_server.
        The get_command_list and execute functions should be implemented in
        each concrete implementation of IWorker. When an IWorker is linked
        to a ServerDaemon, its internal logger function is replaced, but by
        default, the logging function is simply a print to standard out.
"""
# endregion
class IWorker(object):
    def __init__(self):
        self.logger = self.__print
        self.name = None

    # region Method Description
    """
    Method: set_logger
        Description:
            This method allows the ServerDaemon to set the logging method
            of each of its workers as necessary.
        Arguments:
            logging_method: a pointer a method that accepts one argument,
                a string to be logged.
    """
    # endregion
    def set_logger(self, logging_method):
        self.logger = logging_method

    # region Method Description
    """
    Method: __print
        Description:
            An implementation of just a print statement so that the logger
            has a default.
    """
    # endregion
    def __print(self, statement):
        print statement

    # region Method Description
    """
    Method: get_command_list
        Description:
            Method used to return a list of commands that this instance of
            IWorker is responsible for. When the ACC issues commands from this
            list, then, the ServerDaemon will be able to determine that
            this is the worker to be used.
        Returns:
            This method should return a simple list of commands.
            i.e. ['BRICKCAL', 'BRICKMOVE']
    """
    # endregion
    def get_command_list(self):
        raise NotImplementedError

    # region Method Description
    """
    Method: execute
        Description:
            This method will be passed the command from the ACC given that
            the command list related to this IWorker contains the ACC command.
            The implementation of this method should be able to handle every
            command from get_command_list.
        Arguments:
            acc_command: array of command and parameters from the ACC.
    """
    # endregion
    def execute(self, acc_command):
        raise NotImplementedError

    # region Method Description
    """
    Method: stateframe_query
        Description:
            This method is called on by the server to poll data from each
            of the workers. Depending on the worker, the return format
            may be different.
        Returns:
            Generally should return a dictionary with the polled data.
    """
    def stateframe_query(self):
        raise NotImplementedError
