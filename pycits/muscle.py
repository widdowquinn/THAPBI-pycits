#!/usr/bin/env python
#
# Tools for working with MUSCLE
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys

from subprocess import call
from .tools import is_exe


class Muscle(object):
    """Class for working with MUSCLE"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No MUSCLE available at %s (exiting)" %
                               exe_path)
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, seqdir):
        """Run MUSCLE on the passed file"""
        self._outdirname = os.path.join(seqdir, "MUSCLE-aligned_OTUs")
        if not os.path.exists(self._outdirname):
            os.makedirs(self._outdirname)
        infiles = [f for f in os.listdir(seqdir) if
                   os.path.splitext(f)[-1] == '.fasta']
        for fname in infiles:
            self.__build_cmd(seqdir, fname)
            msg = ["Running...", "\t%s" % self._cmd]
            for m in msg:
                self._logger.info(m)
            retcode = call(self._cmd, shell=True)
            if retcode < 0:
                self._logger.error("MUSCLE alignment terminated by " +
                                   "signal %s" % -retcode)
                sys.exit(1)
        return self._outdirname

    def __build_cmd(self, seqdir, fname):
        """Build a command-line for MUSCLE"""
        cmd = ["muscle",
               "-in", os.path.join(seqdir, fname),
               "-out", os.path.join(self._outdirname, fname)]
        self._cmd = ' '.join(cmd)
