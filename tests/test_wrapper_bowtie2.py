#!/usr/bin/env python

"""Tests pycits wrapper for bowtie2.

Bowtie2 does not have native support for gzipped files, so we use
FASTA/FASTQ.
"""

import hashlib
import os
import shutil

from pycits import bowtie2
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal


import subprocess

# INPUT DATA
INDIR = os.path.join("tests", "test_data", "bowtie2")
OUTDIR = os.path.join("tests", "test_out", "bowtie2")
TARGETDIR = os.path.join("tests", "test_targets", "bowtie2")
INFILE = os.path.join(INDIR, "bowtie_test.assembled.fasta")
FA_INDEX = os.path.join(INDIR, "fasta_index")
OUTFILE = os.path.join(OUTDIR, "bowtie2_test.sam")
TARGETFILE = os.path.join(TARGETDIR, "bowtie2_target_mapping.sam")

THREADS = 1


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_bowtie2path():
    """bowtie2 executable is in $PATH"""
    bowtie2.Bowtie2("bowtie2")


def test_bowtie2_cmd():
    """bowtie2 returns correct form of cmd-line"""
    bt2_map = bowtie2.Bowtie2("bowtie2")
    target = ' '.join([str(e) for e in ["bowtie2",
                                        "-x", FA_INDEX,
                                        "-f", INFILE,
                                        "-S", OUTFILE,
                                        "-p", THREADS,
                                        "--very-sensitive",
                                        "--no-unal"]])
    result = bt2_map.run(INFILE, FA_INDEX, OUTFILE, THREADS, fasta=True, dry_run=True)
    assert_equal(result.command, target)


def test_bowtie2_exec_notexist():
    """error thrown if bowtie2 executable does not exist"""
    try:
        obj = bowtie2.Bowtie2(os.path.join(".", "bowtie2"))
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_notexec():
    """error thrown if bowtie2 not executable"""
    try:
        obj = bowtie2.Bowtie2("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_build_exec():
    """bowtie2 maps reads correctly with FASTA test data"""
    bt2_map = bowtie2.Bowtie2("bowtie2")
    result = bt2_map.run(INFILE, FA_INDEX, OUTFILE, THREADS, fasta=True)
    # Test for equality of output and target MD5 files
    # To establish equality, we need to ignore the @PG line that describes the
    # path to bowtie2, as this may differ between systems
    with open(os.path.join(OUTFILE), "r") as outfh:
        with open(TARGETFILE, "r") as tgtfh:
            tgtdata = '\n'.join([l for l in tgtfh.readlines() if
                                 not l.startswith('@PG')])
            outdata = '\n'.join([l for l in outfh.readlines() if
                                 not l.startswith('@PG')])
            assert_equal(tgtdata, outdata)
