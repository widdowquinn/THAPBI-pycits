#!/usr/bin/env python

"""Tests of program to convert fq to fa
This is located in the bin directory. This has been left as a stand alone
tool as it could be useful for other purposes. """

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits import clean_up


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_bio_fq_fa")
TEST_FILE = "pear_test.assembled.fastq.gz"
TEST_OUT = "converted_fq_fa.fasta"
# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "convert_fq_fa",
                      "pear_tests_FASTA_convert.fasta.gz")


def test_convert():
    """convert instantiates with cmd-line if convert_fq_to_fa.py
    is in relative $PATH"""
    cmd = ["python",
           os.path.join("./bin", "convert_fq_to_fa.py"),
           "-h"]
    cmd = ' '.join(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)


def test_convert_exec_notexist():
    """Error thrown if convert_fq_to_fa.py executable does not exist"""
    try:
        # Using the clean up module to test this instead of making a module
        # just for this purpose.
        obj = clean_up.Clean_up(os.path.join("./bin", "convert_fq_to_fa.py"))
    except NotExecutableError:
        return True
    else:
        return False


def test_convert_notexec():
    """Error thrown if convert exe not executable"""
    try:
        obj = clean_up.Clean_up(os.path.join(".", "LICENSE"))
    except NotExecutableError:
        return True
    else:
        return False


def test_convert_exec():
    """Run convert_fq_to_fa.py on test data and compare output to
    precomputed target
    """
    cmd = ["python",
           (os.path.join("./bin", "convert_fq_to_fa.py")),
           "-i",
           (os.path.join(INDIR, TEST_FILE)),
           "-o",
           (os.path.join(OUTDIR, TEST_OUT))]
    cmd = ' '.join(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    result_file = os.path.join(OUTDIR, TEST_OUT)
    with gzip.open(TARGET, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
