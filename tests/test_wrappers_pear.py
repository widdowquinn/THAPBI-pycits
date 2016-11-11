#!/usr/bin/env python

"""Tests of wrapper code in pycits. Trimmomatic"""

import os
import shutil

from pycits import pear
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_pear")

# this will depend on the binarie name of pear
# WE will need to change this!!


def test_pear():
    """pear instantiates with cmd-line if pear is in $PATH"""
    assemble = pear.Pear("pear")


def test_pear_cmd():
    """pear instantiates and returns correct form of cmd-line"""
    assembl = pear.Pear("pear")
    prefix = "DNAMIX_S95_L001"
    outdirname = os.path.join("tests/test_out_pear", "%s_PEAR"
                              % (prefix))
    command = ["pear",
               "-f",
               "tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
               "-r",
               "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
               "--threads",
               "4",
               "-o",
               outdirname]
    target = ' '.join(command)
    assert_equal(assembl.run("tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
                             "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
                             4, "tests/test_out_pear"), target)


def test_pear_exec_notexist():
    """Error thrown if pear executable does not exist"""
    try:
        assemble = pear.Pear(os.path.join(".", "pear"))
    except NotExecutableError:
        return True
    else:
        return False


def test_pear_notexec():
    """Error thrown if pear exe not executable"""
    try:
        assemble = pear.Pear("LICENSE")
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


def test_pear_exec():
    """Run pear on test data"""
    assemble = pear.Pear("pear")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    assemble.run("tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
                 "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
                 4, "tests/test_out_pear/")
    target = os.path.join("tests", "test_targets",
                          "pear_test.assembled.fastq")
    result = os.path.join("tests", "test_out_pear",
                          "DNAMIX_S95_L001_PEAR.assembled.fastq")
    with open(target, "r") as target_fh:
        with open(result, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
