#!/usr/bin/env python
#
# Tools for working with FastQC
#
# FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

import os

from tools import is_exe

class FastQC(object):
    """Class for working with FastQC"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            msg = ["FastQC executable not valid",
                   "QC of reads will not be run"]
            for m in msg:
                self._logger.warning(m)
            self._no_run = True
        self._exe_path = exe_path
            
