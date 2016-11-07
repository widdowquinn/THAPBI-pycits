#!/usr/bin/env python
#
# error_correction * (correct Illumina specific errors)
# http://bioinf.spbau.ru/spades/bayeshammer
# This comes bunduled with SPAdes.
# http://bioinf.spbau.ru/en/content/spades-download-0
# tested with 3.9.0
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from .tools import is_exe


class Error_correction(object):
    """Class for working with Bayes hammer"""

    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            if logger:
                self._logger.error("""No SPADES program in PATH (exiting)
        please download from :
        http://bioinf.spbau.ru/en/content/spades-download-0
        and add to your PATH, or give full path to ...
        """)
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, exe_path, L_reads, R_reads, outdir, threads,
            logger=False):
        """Run SPAdes on the read files"""
        assert L_reads != R_reads, """Oh no,
        I am trying to correct two files that are the same!
        Something has gone wrong in determining the left and right
        read files."""
        self.__build_cmd(exe_path, L_reads, R_reads,
                         str(threads), outdir)
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
                self._logger.error("SPAdes terminated by " +
                                   "signal %s" % pipe.returncode)
            sys.exit(1)
        else:
            if logger:
                self._logger.info("SPAdes.py returned %s"
                                  % pipe.returncode)
        # print ("\n\npipe.args = ", pipe.args, "\n\n")
        return pipe.args

    def __build_cmd(self, exe_path, L_reads, R_reads,
                    threads, outdir):
        """Build a command-line for SPAdes to error correct the
         paired end reads. """
        self._outdirname = os.path.join(outdir)
        cmd = ["python",
               exe_path,
               "-1", L_reads,
               "-2", R_reads,
               "--only-error-correction",
               "--threads", str(threads),
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)
