#!/usr/bin/env python
#
#
# Swarm (clustering)
# https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

# WARNING: this module is not yet tested AT ALL.

import os
import sys

import subprocess
from .tools import is_exe


class Swarm(object):
    """Class for working with SWARM"""
    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        if logger:
            self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            msg = ("""No SWARM program in PATH.""
            please download/ install from:
            https://github.com/torognes/swarm
            and add to your PATH(exiting)""")
            if logger:
                self._logger.warning(msg)
            sys.exit(1)
        self._exe_path = exe_path

    def __build_cmd(self, infname, threads, clustering_threshold,
                    outdir):
        """Build a command-line for swarm"""
        clustering_threshold = str(clustering_threshold)
        self._outdirname = os.path.join(outdir,
                         "swarm_clustering_d%s" %
                         (clustering_threshold))
        cmd = ["swarm",
               "-t", str(threads),
               "-d", clustering_threshold,
               "-o",
               os.path.join(self._outdirname,
                            "swarm_clustering_d%s" %
                            (clustering_threshold)),
               infname] 
        self._cmd = ' '.join(cmd)

    def run(self, infname, threads, clustering_threshold,
                    outdir, logger=False):
        """Run SWARM on the passed file"""
        self.__build_cmd(infname, threads, \
                         clustering_threshold,
                         outdir)
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
                self._logger.error("Swarm terminated by " +
                                   "signal %s" % -retcode)
            sys.exit(1)
        if pipe.returncode == 0:
            if logger:
                self._logger.info(pipe)
        print ("pipe.args = ", pipe.args)
        return (self._outdirname, pipe.args)

