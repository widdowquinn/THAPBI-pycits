#!/usr/bin/env python

"""Test program for a function in the module tools: convert_fq_to_fa
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import convert_fq_to_fa


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


def test_convert_exec():
    """Run convert_fq_to_fa.py on test data and compare output to
    precomputed target
    """
    result_file = os.path.join(OUTDIR, TEST_OUT)
    convert_fq_to_fa((os.path.join(INDIR, TEST_FILE)), result_file)

    with gzip.open(TARGET, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
