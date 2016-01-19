#!/usr/bin/env python
#
# Tools for working with BLAST

import os
import sys

from tools import is_exe

class Blastclust(object):
    """Class for working with blastclust"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No blastclust at %s (exiting)" % exe_path)
            sys.exit(1)
        self._exe_path = exe_path
