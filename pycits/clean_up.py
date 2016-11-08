#!/usr/bin/env python
#
#   module to clean up given files
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys
from .tools import is_exe


from subprocess import call

class Clean_up(object):
    """Class for cleaning up files."""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("error in command")
            sys.exit(1)
        self._exe_path = exe_path


    def run(self, infile):
        """Run clean up on unwanted files"""
        self.__build_cmd(infile)
        msg = ["Running...", "\t%s" % self._cmd]
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("clean up terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("clean_up.py returned %s" % retcode)

    def __build_cmd(self, infile):
        """Build a command-line to clean up files """
        cmd = ["rm",
               infile]
        self._cmd = ' '.join(cmd)
