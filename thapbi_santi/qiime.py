#!/usr/bin/env python
#
# Tools for working with QIIME scripts
#
# QIIME: http://qiime.org/

import os
import sys

from subprocess import call
from tools import is_exe

class Pick_Otus(object):
    """Class for working with pick_otus.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No pick_otus.py script (exiting)")
            sys.exit(1)
        self._exe_path = exe_path

class Pick_Closed_Ref_Otus(object):
    """Class for working with pick_closed_reference_otus.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No pick_closed_reference_otus.py script " +
                               "(exiting)")
            sys.exit(1)
        self._exe_path = exe_path
            
class Join_Paired_Ends(object):
    """Class for working with join_paired_ends.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No join_paired_ends.py script (exiting)")
            sys.exit(1)
        self._exe_path = exe_path
            
    def run(self, infnames, outdir):
        """Run joined_paired_ends.py on the passed file"""
        self.__build_cmd(infnames, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("join_paired_ends.py terminated by signal %s" %
                               -retcode)
            sys.exit(1)
        else:
            self._logger.info("join_paired_ends.py returned %s" % retcode)
        return self._outfname

    def __build_cmd(self, infnames, outdir):
        """Build a command-line for join_paired_ends.py"""
        f1, f2 = tuple(infnames)
        cmd = ["join_paired_ends.py",
               "-f", f1,
               "-r", f2,
               "-o", outdir]
        self._cmd = ' '.join(cmd)
        self._outfname = os.path.join(outdir, "fastqjoin.join.fastq")
