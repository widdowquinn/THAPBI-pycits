#!/usr/bin/env python

"""Tests pycits wrapper for muscle."""

import os
import shutil

from pycits import muscle
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_muscle")
FASTA = os.path.join(INDIR, "dedup_test.fasta")

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "muscle",
                      "muscle_tests.fasta")


def test_muscle():
    """muscle instantiates with cmd-line if muscle is in $PATH"""
    muscle.muscle("muscle")


def test_muscle_cmd():
    """muscle instantiates, runs and returns correct form of cmd-line"""
    obj = muscle.Muscle("muscle")
    target = ' '.join(["muscle", FASTA, "-o", OUTDIR])
    qc = obj.run(FASTA, OUTDIR)
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


# there isnt any reason to test this ...
def test_muscle_exec():
    """Run muscle on test data and compare output to precomputed target"""
    obj = muscle.Muscle("muscle")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = obj.run(READS1, OUTDIR)
    with open(TARGET, "r") as target_fh:
        with open(result.muscle_html, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
