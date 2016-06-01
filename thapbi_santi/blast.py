#!/usr/bin/env python
#
# Tools for working with BLAST
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys

from subprocess import Popen, PIPE
from tools import is_exe


class Blastclust(object):
    """Class for working with blastclust"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No blastclust at %s (exiting)" % exe_path)
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, outdir, threads):
        """Run blastclust on the passed file"""
        self.__build_cmd(infnames, outdir, threads)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        pipe = Popen(self._cmd, shell=True, stdout=PIPE)
        if pipe.wait() != 0:
            self._logger.error("blastclust generated some errors")
            sys.exit(1)
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
