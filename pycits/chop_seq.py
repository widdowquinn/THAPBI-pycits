#!/usr/bin/env python
#
# Tools for chopping the assembled seq
# and renaming the databse name with the
# abundancevalue. We assume THERE ARE 1
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from .tools import is_exe


class Chop_seq(object):
    """Class for working with fasta seq
    chopper script"""
    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        if logger:
            self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("""No
            trim_fasta_file.py
            available at %s (exiting)""" %
                               exe_path)
            if logger:
                self._logger.warning(msg)
            self._no_run = True
        self._exe_path = exe_path

    def __build_cmd(self, exe_path, fasta, database,
                    left, right, fasta_out):
        """Build a command-line for
        trim_fasta_file.py"""
        cmd = ["python",
               exe_path,
               "-f", fasta,
               "--database", database,
               "-l", str(left),
               "-r", str(right),
               "-o", fasta_out]
        self._cmd = ' '.join(cmd)

    def run(self, exe_path, fasta, database, left, right,
            fasta_out, logger=False):
        """function to run the command
        trim_fasta_file.py script"""
        self.__build_cmd(exe_path, fasta, database, left,
                         right, fasta_out, barcode)
        msg = ["Running...", "\t%s" % self._cmd]
        # for m in msg:
        #    self._logger.info(m)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if pipe.returncode != 0:
            if logger:
                self._logger.error("""
                trim_fasta_file.py
                generated some errors""")
            sys.exit(1)
        if pipe.returncode == 0:
            if logger:
                self._logger.info(pipe)
        print ("pipe.args = ", pipe.args)
        return (pipe.args)
