#!/usr/bin/env python
#
#
# Swarm (clustering)
# https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

# WARNING: this module is not yet tested AT ALL.

import os
import sys

import subprocess
from .tools import is_exe, NotExecutableError


class SwarmError(Exception):
    """Exception raised when swarm fails"""
    def __init__(self, message):
        self.message = message


class Swarm(object):
    """Class for working with SWARM"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path


    def run(self, infname, outdir, threads, threshold, dry_run=False):
        """Run swarm to cluster sequences in the passed file

        - infname    - path to sequences for clustering
        - outdir     - output directory for clustered output
        - threads    - number of threads for swarm to use
        - threshold  - clustering threshold for swarm (-d option)
        - dry_run    - if True returns cmd-line but does not run
        """
        self.__build_cmd(infname, threads, threshold, outdir)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True) 
        return (self._outfname, pipe.stdout.decode('utf-8'))


    def __build_cmd(self, infname, threads, threshold, outdir):
        """Build a command-line for swarm"""
        self._outfname = os.path.join(outdir, "swarm.out")
        cmd = ["swarm",
               "-t", str(threads),
               "-d", str(threshold),
               "-o", self._outfname,
               infname] 
        self._cmd = ' '.join(cmd)


