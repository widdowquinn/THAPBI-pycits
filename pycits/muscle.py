#!/usr/bin/env python
#
# Tools for working with MUSCLE
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Flash class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# outfile - this is the out file
# stderr
Results = namedtuple("Results", "command outfile " +
                     "stdout stderr")

class MuscleError(Exception):
    """Exception raised when Muscle fails"""
    def __init__(self, message):
        self.message = message

class Muscle(object):
    """Class for working with MUSCLE"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path


    def run(self, seqdir, dry_run=False):
        """Run MUSCLE on the passed folder"""
        self._outdirname = os.path.join(seqdir, "MUSCLE-aligned_OTUs")
        if not os.path.exists(self._outdirname):
            os.makedirs(self._outdirname)
        infiles = [f for f in os.listdir(seqdir) if
                   os.path.splitext(f)[-1] == '.fasta']
        for fname in infiles:
            self.__build_cmd(seqdir, fname)
            pipe = subprocess.run(self._cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, seqdir, fname):
        """Build a command-line for MUSCLE"""
        self._outfnames = os.path.join(self._outdirname, fname)
        cmd = ["muscle",
               "-in", os.path.join(seqdir, fname),
               "-out", os.path.join(self._outdirname, fname)]
        self._cmd = ' '.join(cmd)
        
