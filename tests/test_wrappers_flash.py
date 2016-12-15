#!/usr/bin/env python

"""Tests of wrapper code in pycits: Flash"""

import os
import shutil

from pycits import flash
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_flash")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = os.path.split(READS1)[-1].split("_R")[0]

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "flash",
                      "flash_test.assembled.fastq")

def test_flash():
    """flash instantiates with cmd-line if flash is in $PATH"""
    assemble = flash.Flash("flash")

def test_flash_cmd():
    """flash instantiates, runs and returns correct form of cmd-line"""
    obj = flash.Flash("flash")
    target = ' '.join(["flash", "--max-overlap", "250",
                       "--threads", "4", "-o", os.path.join(OUTDIR, PREFIX),
                       READS1, READS2])
    assert_equal(obj.run(READS1, READS2, 4, OUTDIR, PREFIX, dry_run=True),
                 target)

def test_flash_exec_notexist():
    """Error thrown if flash executable does not exist"""
    try:
        obj = flash.Flash(os.path.join(".", "flash"))
    except NotExecutableError:
        return True
    else:
        return False

def test_flash_notexec():
    """Error thrown if flash exe not executable"""
    try:
        obj = flash.Flash("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_flash_exec():
    """Run flash on test data and compare output to precomputed target
    this often does not produce the same result twice. So we may have
    to disable this test"""
    obj = flash.Flash("flash")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = obj.run(READS1, READS2, 4, OUTDIR, PREFIX, dry_run=False)
    with open(TARGET, "r") as target_fh:
        with open(result.outfileassembled, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
