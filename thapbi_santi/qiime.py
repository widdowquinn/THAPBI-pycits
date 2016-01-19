#!/usr/bin/env python
#
# Tools for working with QIIME scripts
#
# QIIME: http://qiime.org/

import os
import sys

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
            
