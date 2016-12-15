#!/usr/bin/env python
#
# FLASH * (assemble overlapping reads)
# https://ccb.jhu.edu/software/FLASH/
# follow this link to get the binaries.
# https://sourceforge.net/projects/flashpage/
# ?source=typ_redirect
# This has been tested with FLASH-1.2.11

#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess

from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Pear class returned values
Results = namedtuple("Results", "command outfileassembled outfilediscarded " +
                     "outfileunassmbledfwd outfileunassembledrev " +
                     "stdout stderr")


class FlashError(Exception):
    """Exception raised when flash fails"""
    def __init__(self, message):
        self.message = message

class Flash(object):
    """class for working with paired end read assembly tool Flash"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, lreads, rreads, threads, outdir, prefix, dry_run=False):
        """Run Flash to merge passed read files

        - lreads    - forward reads
        - rreads    - reverse reads
        - threads   - number of threads for pear to use
        - outdir    - output directory for merged output
        - prefix    - file prefix for flash output
        - logger    - stream to write messages
        - dry_run   - if True, returns cmd-line but does not run

        Returns a tuple of output filenames, and the STOUT returned by the
        flash run.
        """
        assert(lreads != rreads)  # We don't want to merge the same file
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
