#!/usr/bin/env python

import os

from subprocess import check_output

def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        exefile = check_output(["which", filename]).strip()
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)
