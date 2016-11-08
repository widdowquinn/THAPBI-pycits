#!/usr/bin/env python
#
# Tools for deduplicating reads
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from .tools import is_exe

class Deduplicate(object):
    """Class for working with deduplicate script"""
    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        if logger:
            self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("""No deduplicate_rename.py
                               available at %s (exiting)""" %
                               exe_path)
            if logger:
                self._logger.warning(msg)
            self._no_run = True
        self._exe_path = exe_path
        
    def __build_cmd(self, exe_path, fasta, database,
                    out, logger=False):
        """Build a command-line for deduplicate_rename.py"""
        cmd = ["python",
               exe_path,
               "--fasta", fasta,
               "--database", database,
               "-o", out]
        self._cmd = ' '.join(cmd)

    def run(self, exe_path, fasta, database,
                    out, logger=False):
        """deduplicate fasta file"""
        self.__build_cmd(exe_path, fasta, database,
                         out, logger)
        msg = ["Running...", "\t%s" % self._cmd]
        # for m in msg:
        #    self._logger.info(m)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)            
        if pipe.returncode != 0:
            if logger:
                self._logger.error("""deduplicate_rename
                                   generated some errors""")
            sys.exit(1)
        if pipe.returncode == 0:
            if logger:
                self._logger.info(pipe)
        print ("pipe.args = ", pipe.args)
        return (self._outdirname, pipe.args)




