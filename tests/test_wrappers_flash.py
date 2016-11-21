#!/usr/bin/env python

"""Tests of wrapper code in pycits. Trimmomatic"""

import os
import shutil

from pycits import flash
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_flash")

# this will depend on the binarie name of flash
# WE will need to change this!!


def test_flash():
    """flash instantiates with cmd-line if flash is in $PATH"""
    assemble = flash.Flash("flash")


def test_flash_cmd():
    """flash instantiates and returns correct form of cmd-line"""
    assembl = flash.Flash("flash")
    prefix = "DNAMIX_S95_L001"
    outdirname = os.path.join("tests/test_out_flash", "%s_flash"
                              % (prefix))
    command = ["flash",
               "--max-overlap", "250",
               "--threads",
               "4",
               "-o",
               outdirname,
               "tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
               "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",]
    target = ' '.join(command)
    assert_equal(assembl.run("tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
                             "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
                             4, "tests/test_out_flash"), target)


def test_flash_exec_notexist():
    """Error thrown if flash executable does not exist"""
    try:
        assemble = flash.Flash(os.path.join(".", "flash"))
    except NotExecutableError:
        return True
    else:
        return False


def test_flash_notexec():
    """Error thrown if flash exe not executable"""
    try:
        assemble = flash.Flash("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def build_diff_cmd(infname1, infname2):
    """Build a command-line for diff"""
    cmd = ["diff",
           infname1,
           infname2]
    build_diff_cmd = ' '.join(cmd)
    return build_diff_cmd

# cannot get the same result twice. 
##def test_flash_exec():
##    """Run flash on test data"""
##    assemble = flash.Flash("flash")
##    try:
##        shutil.rmtree(OUTDIR)
##    except FileNotFoundError:
##        pass
##    os.makedirs(OUTDIR, exist_ok=True)
##    assemble.run("tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
##                 "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
##                 4, "tests/test_out_flash/")
##    target = os.path.join("tests", "test_targets",
##                          "DNAMIX_S95_L001_flash.extendedFrags.fastq")
##    result = os.path.join("tests", "test_out_flash",
##                          "DNAMIX_S95_L001_flash.extendedFrags.fastq")
##    with open(target, "r") as target_fh:
##        with open(result, "r") as test_fh:
##            assert_equal(target_fh.read(), test_fh.read())
