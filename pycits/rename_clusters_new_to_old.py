#!/usr/bin/env python
#
# Tools for renaming the swarm clustered output
# with the old names, replacing the "new names"
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from .tools import is_exe


class Rename(object):
    """Class for working with rename script"""
    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        if logger:
            self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("""No
            parse_clusters_new_to_old_name.py
            available at %s (exiting)""" %
                               exe_path)
            if logger:
                self._logger.warning(msg)
            self._no_run = True
        self._exe_path = exe_path

    def __build_cmd(self, exe_path, infile, ITS_database,
                    database,
                    out, logger=False):
        """Build a command-line for
        parse_clusters_new_to_old_name.py"""
        cmd = ["python",
               exe_path,
               "-i", infile,
               "--seq_db", ITS_database,
               "--database", database,
               "-o", out]
        self._cmd = ' '.join(cmd)

    def run(self, exe_path, infile, ITS_database, database,
            out, logger=False):
        """deduplicate fasta file"""
        self.__build_cmd(exe_path, infile, ITS_database,
                         database,
                         out, logger)
        msg = ["Running...", "\t%s" % self._cmd]
        # for m in msg:
        #    self._logger.info(m)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if pipe.returncode != 0:
            if logger:
                self._logger.error("""
                parse_clusters_new_to_old_name.py
                generated some errors""")
            sys.exit(1)
        if pipe.returncode == 0:
            if logger:
                self._logger.info(pipe)
        print ("pipe.args = ", pipe.args)
        return (pipe.args)
