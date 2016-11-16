#!/usr/bin/env python
#
# Tools for working with BLAST
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys

from subprocess import Popen, PIPE
from .tools import is_exe, NotExecutableError


class BlastclustError(Exception):
    """Exception raised when blastclust fails"""
    def __init__(self, message):
        self.message = message


class Blastclust(object):
    """Class for working with blastclust"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infnames, outdir, threads, dry_run=False):
        """Run blastclust to cluster sequences in the passed file.

        - infnames - input filenames for clustering
        - outdir   - output directory for clustered output
        - threads  - number of threads for blastclust to use
        - dry_run  - if True, returns cmd-line but does not run
        """
        self.__build_cmd(infnames, outdir, threads)
        if dry_run:
            return self._cmd
        pipe = Popen(self._cmd, shell=True, stdout=PIPE)
        if pipe.wait() != 0:
            msg = "blastclust generated some errors"
            raise BlastclustError(msg)
        return (self._outfname, pipe.stdout.readlines())

    def __build_cmd(self, infname, outdir, threads):
        """Build a command-line for blastclust"""
        self._outfname = os.path.join(outdir,
                                      os.path.split(infname)[-1] +
                                      ".blastclust99.lst")
        cmd = ["blastclust",
               "-L", "0.90", "-S", "99", "-a", str(threads), "-p", "F",
               "-i", infname,
               "-o", self._outfname]
        self._cmd = ' '.join(cmd)
