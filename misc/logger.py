#!/usr/bin/env python

'''
@author: Dr. Abdullah H. Sahyoun, TRON gGmbH
@author: Patrick Sorn, TRON gGmbH

Copyright (c) 2016 TRON gGmbH

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''


'''This module is a wrapper around the built-in logging module.'''

import os
import time

DEBUG = 20
INFO = 30
ERROR = 40

class Logger(object):
    def __init__(self, infile):
        '''This function initializes the logger with log level INFO.'''
        self.logfile = infile
        self.level = 30
        self.levelString = ""

    def setLevel(self, level):
        '''This method sets the level of logging. There are three types of log levels available: DEBUG, INFO, ERROR.'''
        self.level = level
        if level == DEBUG:
            self.levelString = "DEBUG"
        elif level == INFO:
            self.levelString = "INFO"
        elif level == ERROR:
            self.levelString = "ERROR"

    def log(self, level, message):
        '''This method logs a specific message on the specified log level.'''
        log = None
        log_line = time.strftime("%Y-%m-%d %H:%M:%S ") + self.levelString + " " + str(message)
        if os.path.exists(self.logfile):
            log = open(self.logfile,"a")
        else:
            log = open(self.logfile,"w")
        print(log_line)
        log.write(log_line+"\n")
        log.close()

    def debug(self, message):
        '''This method logs a specific message on the DEBUG log level.'''
        if self.level > 20:
            return
        log_line = time.strftime("%Y-%m-%d %H:%M:%S ") + "DEBUG " + str(message)
        if os.path.exists(self.logfile):
            log = open(self.logfile,"a")
        else:
            log = open(self.logfile,"w")
        print(log_line)
        log.write(log_line+"\n")
        log.close()

    def info(self, message):
        '''This method logs a specific message on the INFO log level.'''
        if self.level > 30:
            return
        log_line = time.strftime("%Y-%m-%d %H:%M:%S ") + "INFO " + str(message)
        if os.path.exists(self.logfile):
            log = open(self.logfile,"a")
        else:
            log = open(self.logfile,"w")
        print(log_line)
        log.write(log_line+"\n")
        log.close()

    def error(self, message):
        '''This method logs a specific message on the ERROR log level.'''
        log_line = time.strftime("%Y-%m-%d %H:%M:%S ") + "ERROR " + str(message)
        if os.path.exists(self.logfile):
            log = open(self.logfile,"a")
        else:
            log = open(self.logfile,"w")
        print(log_line)
        log.write(log_line+"\n")
        log.close()


def main():
    logger = Logger('test.log')
    logger.log('TEST','test123')
    logger.debug('error occured')


if __name__ == '__main__':
    main()
