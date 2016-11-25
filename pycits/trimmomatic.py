#!/usr/bin/env python
#
# Tools for working with trimmomatic
#
# trimmomatic:
# http://www.bioinformatics.babraham.ac.uk/projects/trimmomatic/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe


import os
import subprocess
import sys

from collections import namedtuple

from .tools import is_exe, NotExecutableError


# factory class for Trimmomatic class returned values
Results = namedtuple("Results", "command " +
                     "outfileR1paired outfileR1unpaired " +
                     "outfileR2paired outfileR2unpaired " +
                     "stdout stderr")

# factory class for Trimmomatic parameters
Parameters = namedtuple("Parameters", "threads")
Parameters.__new__.__defaults__ = (1,)

# factory class for Trimmomatic trimming steps
Steps = namedtuple("Steps", "ILLUMINACLIP LEADING HEADCROP TRAILING MINLEN " +
                   "SLIDINGWINDOW")
Steps.__new__.__defaults__ = (None, 3, 0, 3, 70, "4:25")


class TrimmomaticError(Exception):
    """Exception raised when Trimmomatic fails"""
    def __init__(self, message):
        self.message = message


class Trimmomatic(object):
    """Class for working with trimmomatic"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, lreads, rreads, outdir, prefix, phred, parameters, steps,
            dry_run=False):
        """Run trimmomatic to trim reads in the passed files, on quality

        - lreads      - forward reads
        - rreads      - reverse reads
        - outdir      - directory to write trimmed output
        - prefix      - prefix string for output files
        - phred       - phred33|phred64
        - parameters  - namedtuple for Trimmomatic parameters
        - steps       - namedtuple for Trimmomatic steps
        - dry_run     - returns only the command to be run, if True

        Returns namedtuple with form
          "command outfileR1paired outfileR1unpaired outfileR2paired
           outfileR2unpaired stdout stderr"

        TODO: accept arbitrary options provided by the user
        """
        assert(lreads != rreads)
        self.__build_cmd(lreads, rreads, outdir, prefix, phred,
                         parameters, steps)
        if dry_run:
            return self._cmd
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, *self._outfnames,
                          pipe.stdout, pipe.stderr)
        return results

    def __build_cmd(self, lreads, rreads, outdir, prefix, phred,
                    parameters, steps):
        """Build a command-line for trimmomatic

        - lreads      - forward reads
        - rreads      - reverse reads
        - outdir      - directory to write trimmed output
        - prefix      - prefix string for output files
        - phred       - phred33|phred64
        - parameters  - namedtuple for Trimmomatic parameters
        - steps       - namedtuple for Trimmomatic steps
        - dry_run     - returns only the command to be run, if True
        """
        self._outfnames = [os.path.join(outdir, prefix + suffix) for suffix in
                           ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
                            "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]
        params = ["-{0} {1}".format(k, v) for (k, v) in
                  parameters._asdict().items() if v is not None]
        steps = ["{0}:{1}".format(k, v) for (k, v) in
                 steps._asdict().items() if v is not None]
        cmd = ["trimmomatic", "PE", "-{0}".format(phred), *params,
               lreads, rreads, *self._outfnames, *steps]
        self._cmd = ' '.join(cmd)
