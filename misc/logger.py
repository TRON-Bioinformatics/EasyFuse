#!/usr/bin/env python

"""
Process logging
Either creates (if not existing) or appends message lines to a file
Message lines are prefixed with a time stamp and keywords (logging level) for the
type of info presented in the line
For convenience, the individual log levels can be accessed via their own methods

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

import time

NONE = 10
DEBUG = 20
INFO = 30
ERROR = 40

class Logger(object):
    """Logging of processing steps"""
    def __init__(self, infile):
        """Variable initialization"""
        self._file = infile
        self.logfile = infile
        self.level = INFO
        self.level_string = ""

    def _set_level(self, level):
        """Set lvl string at the start of logging line"""
        self.level = level
        if level == NONE:
            self.level_string = "NONE"
        elif level == DEBUG:
            self.level_string = "DEBUG"
        elif level == INFO:
            self.level_string = "INFO"
        elif level == ERROR:
            self.level_string = "ERROR"

    # writing message to file
    def _log(self, message):
        """General logging method"""
        log_line = ""
        if self.level_string == "NONE":
            log_line = str(message)
        else:
            log_line = "{0} {1} {2}".format(time.strftime("%Y-%m-%d %H:%M:%S"), self.level_string, str(message)) # pylint: disable=line-too-long
        # a+ either creates a new file or appends to an existing one
        # note - urla: I think direct reading after writing to such a file is not possible
        #              w/o seeking back to the start of the file (not 100% sure).
        #              But this should not be relevant for a logger
        with open(self.logfile, "a+") as log:
            log.write(log_line + "\n")

    def debug(self, message):
        """Logging of debug messages"""
        self._set_level(DEBUG)
        self._log(message)

    def info(self, message):
        """Logging of info messages"""
        self._set_level(INFO)
        self._log(message)

    def error(self, message):
        """Logging of error messages"""
        self._set_level(ERROR)
        self._log(message)

    def no_lvl(self, message):
        """Logging with time/lvl prefix"""
        self._set_level(NONE)
        self._log(message)

    def get_path(self):
        """return the fullpath string of this config file"""
        return self._file
