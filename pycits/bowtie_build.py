#!/usr/bin/env python
#
# Tools for working with Bowtie2 build
#
# http://bowtie-bio.sourceforge.net/bowtie2/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for bowtie build class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# index - this is the index file generated.
# stderr
Results = namedtuple("Results", "command index stdout stderr")


class Bowtie2_Build(Exception):
    """Exception raised when Bowtie2_Build fails"""
    def __init__(self, message):
        self.message = message


class Bowtie2_Build(object):
    """Class for working with bowtie2-build"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infname, index_out, dry_run=False):
        """
        Bowtie2 is btter for longer reads, than bowtie 1.
        Run bowtie2-build on the passed file
        infnames    -  the fasta to index
        index_out   -  the out index prefix
        """
        self.__build_cmd(infname, index_out)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfname, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infname, index_out):
        """Build a command-line for bowtie2-build"""
        self._outfname = index_out
        cmd = ["bowtie2-build",
               "--quiet",
               "-f",
               infname,
               self._outfname]
        self._cmd = ' '.join(cmd)
