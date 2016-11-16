#!/usr/bin/env python

"""Tests of wrapper code in pycits."""

import os
import shutil

from pycits import blast
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal


INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_blastclust")


def test_blastclust():
    """Blastclust instantiates with cmd-line if blastclust is in $PATH"""
    bc = blast.Blastclust("blastclust")


def test_blastclust_cmd():
    """Blastclust instantiates and returns correct form of cmd-line"""
    bc = blast.Blastclust("blastclust")
    target = ' '.join(["blastclust -L 0.90 -S 99 -a 4 -p F",
                       "-i trimmed.fasta",
                       "-o test_out/trimmed.fasta.blastclust99.lst"])
    assert_equal(bc.run("trimmed.fasta", "test_out", 4, dry_run=True),
                 target)


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


def test_blastclust_exec():
    """Run blastclust on test data"""
    bc = blast.Blastclust("blastclust")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = bc.run(os.path.join(INDIR, "trimmed.fasta"), OUTDIR, 4)
    target = os.path.join("tests", "test_targets",
                          "target_trimmed.fasta.blastclust99.lst")
    with open(target, "r") as target_fh:
        with open(result[0], "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
