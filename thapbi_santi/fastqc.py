#!/usr/bin/env python
#
# Tools for working with FastQC
#
# FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

import os
import sys

from subprocess import Popen, PIPE
from tools import is_exe


class FastQC(object):
    """Class for working with FastQC"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            msg = ["FastQC executable not valid",
                   "QC of reads will not be run"]
            for m in msg:
                self._logger.warning(m)
            self._no_run = True
        self._exe_path = exe_path

    def run(self, infnames, outdir):
        """Run fastqc on the passed file"""
        self.__build_cmd(infnames, outdir)
        if not os.path.exists(self._outdirname):
            self._logger.info("Creating output directory: %s" %
                              self._outdirname)
            os.makedirs(self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        pipe = Popen(self._cmd, shell=True, stdout=PIPE)
        if pipe.wait() != 0:
            self._logger.error("fastqc generated some errors")
            sys.exit(1)
        return (self._outdirname, pipe.stdout.readlines())

    def __build_cmd(self, infname, outdir):
        """Build a command-line for fastqc"""
        self._outdirname = os.path.join(outdir, "FastQC_output")
        cmd = ["fastqc",
               "-i", infname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)
