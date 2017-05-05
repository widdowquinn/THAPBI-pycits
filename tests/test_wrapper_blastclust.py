#!/usr/bin/env python

"""Tests pycits wrapper for blastclust.
Also tests the conversion to another format"""

import os
import shutil

from pycits import blast
from pycits.tools import NotExecutableError
from pycits.tools import reformat_swarm_cls

from nose.tools import nottest, assert_equal, with_setup

THREADS = 2

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "blastclust")
OUTDIR = os.path.join("tests", "test_out", "blastclust")

# TARGET OUTPUT
TARGET = os.path.join("tests", "test_targets", "blastclust",
                      "target_trimmed.fasta.blastclust99.lst")

TARGET_FOR_R = os.path.join("tests", "test_targets", "blastclust",
                            "trimmed.fasta.blastclust99.Rout")


# The setup_outdir() function decorates functions by creating the appropriate
# output directory tree
def setup_outdir():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_blastclust():
    """blastclust is in $PATH"""
    blast.Blastclust("blastclust")


def test_blastclust_cmd():
    """Blastclust instantiates and returns correct form of cmd-line"""
    bc = blast.Blastclust("blastclust")
    target = ' '.join(["blastclust", "-L 0.90", "-S 99", "-a 4", "-p F",
                       "-i trimmed.fasta",
                       "-o", os.path.join("tests", "test_out", "blastclust",
                                          "trimmed.fasta.blastclust99.lst")])
    result = bc.run("trimmed.fasta", OUTDIR, 4, dry_run=True)
    assert_equal(result.command, target)


def test_blastclust_exec_notexist():
    """Error thrown if executable does not exist"""
    try:
        bc = blast.Blastclust(os.path.join(".", "blastclust"))
    except NotExecutableError:
        return True
    else:
        return False


def test_blastclust_notexec():
    """Error thrown if blastclust exe not executable"""
    try:
        bc = blast.Blastclust("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


@with_setup(setup_outdir)
def test_blastclust_exec():
    """Run blastclust on test data and compare output to precomputed target"""
    bc = blast.Blastclust("blastclust")
    result = bc.run(os.path.join(INDIR, "trimmed.fasta"), OUTDIR, THREADS)
    with open(TARGET, "r") as target_fh:
        with open(result.outfilename, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


# use the function for converting swarm to "R" format here on this data
def test_convert_exec():
    """run the function reformat_swarm_cls
    to compare the results against
    pregenerated results"""

    result_file_r = os.path.join(OUTDIR,
                                 "trimmed.fasta.blastclust99.Rout")

    # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
    reformat_swarm_cls(os.path.join(OUTDIR,
                                    "trimmed.fasta.blastclust99.lst"),
                       os.path.join(INDIR, "phy_db_forblastclust.fasta"),
                       os.path.join(INDIR, "trimmed.fasta"),
                       result_file_r,
                       False)

    with open(TARGET_FOR_R, "rt") as target_fh:
        data = target_fh.read().split("\n")
        with open(result_file_r, "r") as test_fh:
            result = test_fh.read().split("\n")

            assert_equal(result, data)
