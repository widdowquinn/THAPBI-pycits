#!/usr/bin/env python

"""Tests pycits wrapper for muscle."""

import os
import shutil

from pycits import muscle
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTFOLDER = os.path.join("tests", "test_out_muscle")
OUTFILE = "muscle_test_run.fasta"
FASTA = os.path.join(INDIR, "dedup_test.fasta")

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "muscle",
                      "muscle_tests.fasta")


def test_folder_exists():
    """tests to see if the folder exists for the output folder"""
    if not os.path.exists(OUTFOLDER):
        os.makedirs(OUTFOLDER)


def test_muscle():
    """muscle instantiates with cmd-line if muscle is in $PATH"""
    muscle.Muscle("muscle")


def test_muscle_cmd():
    """muscle instantiates, runs and returns correct form of cmd-line"""
    obj = muscle.Muscle("muscle")
    target = ' '.join(["muscle", "-in", FASTA, "-out",
                       (os.path.join(OUTFOLDER, OUTFILE))])
    qc = obj.run(FASTA, OUTFILE, OUTFOLDER)
    assert_equal(qc.command,
                 target)


def test_muscle_exec_notexist():
    """Error thrown if muscle executable does not exist"""
    try:
        obj = muscle.Muscle(os.path.join(".", "muscle"))
    except NotExecutableError:
        return True
    else:
        return False


def test_muscle_notexec():
    """Error thrown if muscle exe not executable"""
    try:
        obj = muscle.Muscle("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_muscle_exec():
    """Run muscle on test data and compare output to precomputed target"""
    obj = muscle.Muscle("muscle")
    try:
        shutil.rmtree(OUTFILE)
    except FileNotFoundError:
        pass
    os.makedirs(OUTFILE, exist_ok=True)
    result = obj.run(FASTA, OUTFILE, OUTFOLDER)
    with open(TARGET, "r") as target_fh:
        with open(result.muscle_outfile, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
