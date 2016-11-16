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
import sys

import subprocess
from .tools import is_exe, NotExecutableError


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

    def run(self, lreads, rreads, threads, outdir, prefix, adapters,
            dry_run=False):
        """Run trimmomatic to trim reads in the passed files, on quality

        - lreads    - forward reads
        - rreads    - reverse reads
        - threads   - number of threads for Trimmomatic to use
        - outdir    - directory to write trimmed output
        - prefix    - prefix string for output files
        - adapters  - path to adapters to be used with ILLUMINACLIP
        - dry_run   - returns only the command to be run, if True

        TODO: accept arbitrary options provided by the user
        """
        assert(lreads != rreads)
        self.__build_cmd(lreads, rreads, threads, outdir, prefix, adapters)
        if dry_run:
            return self._cmd
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        if pipe.returncode != 0:
            msg = "Trimmomatic returned non-zero: {0}".format(pipe.returncode)
            raise TrimmomaticError('\n'.join(msg, pipe.stdout, pipe.stderr))
        return (self._outfnames, pipe.stdout.decode('utf-8'))

    def __build_cmd(self, lreads, rreads, threads, outdir, prefix, adapters):
        """Build a command-line for trimmomatic"""
        self._outfnames = [os.path.join(outdir, prefix + suffix) for suffix in
                           ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
                            "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]
        cmd = ["trimmomatic", "PE", "-threads {0}".format(threads),
               "-phred33", lreads, rreads,
               *self._outfnames,
               "ILLUMINACLIP:{0}:2:30:10".format(adapters),
               "LEADING:3", "HEADCROP:40", "TRAILING:3",
               "SLIDINGWINDOW:4:25", "MINLEN:70"]
        self._cmd = ' '.join(cmd)
