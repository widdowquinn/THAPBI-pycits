#!/usr/bin/env python
#
# FLASH * (assemble overlapping reads)
# https://ccb.jhu.edu/software/FLASH/
# follow this link to get the binaries.
# https://sourceforge.net/projects/flashpage/
# ?source=typ_redirect

#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from .tools import is_exe


class Flash(object):
    """Class for working with flash"""

    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            if logger:
                self._logger.error("""No flash program in PATH (exiting)
         Download, decompress, type make. Add current directory to you PATH
         . This will create a copy of
        flash in that name. Make sure the PATH to this bin directory is in
        your PATH.
        This has been tested with FLASH-1.2.11

        """)
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, L_reads, R_reads, threads, outdir,
            logger=False):
        """Run flash on the read files"""
        assert L_reads != R_reads, """Oh no,
        I am trying to assemble two files that are the same!
        Something has gone wrong in determining the left and right
        read files."""
        self.__build_cmd(L_reads, R_reads, str(threads), outdir)
        if not os.path.exists(self._outdirname):
            if logger:
                self._logger.info("Creating output directory: %s" %
                                  self._outdirname)
            os.makedirs(self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            if logger:
                self._logger.info(m)

        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if pipe.returncode != 0:
            if logger:
                self._logger.error("flash terminated by " +
                                   "signal %s" % pipe.returncode)
            sys.exit(1)
        else:
            if logger:
                self._logger.info("flash.py returned %s"
                                  % pipe.returncode)
        # print ("\n\npipe.args = ", pipe.args, "\n\n")
        return pipe.args

    def __build_cmd(self, L_reads, R_reads, threads, outdir):
        """Build a command-line for flash to assemble the
        overlapping paired end reads. """
        prefix = L_reads.split("_R")[0]
        prefix = prefix.split("/")[-1]
        self._outdirname = os.path.join(outdir, "%s_flash" % (prefix))
        cmd = ["flash",
               "--max-overlap", "250",
               "--threads", str(threads),
               "-o", self._outdirname,
               L_reads,
               R_reads]
        self._cmd = ' '.join(cmd)
