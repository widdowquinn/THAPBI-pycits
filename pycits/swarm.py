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

from subprocess import call
from .tools import is_exe


class Swarm(object):
    """Class for working with SWARM"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("""No SWARM program in PATH.""
            please download/ install from:
            https://github.com/torognes/swarm
            and add to your PATH(exiting)""")
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
               "--append-abundance", "1",
               "-d", clustering_threshold,
               "-o",
               os.path.join(self._outdirname,
                            "swarm_clustering_d%s" %
                            (clustering_threshold)),
               infname] 
        self._cmd = ' '.join(cmd)

    def run(self, infname, threads, clustering_threshold,
                    outdir):
        """Run SWARM on the passed file"""
        self.__build_cmd(infname, threads, \
                         clustering_threshold,
                         outdir)
        if not os.path.exists(self._outdirname):
            self._logger.info("Creating output directory: %s" %
                              self._outdirname)
            os.makedirs(self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("Swarm terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("SWARM returned %s" % retcode)
        return self._outdirname

