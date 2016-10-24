#!/usr/bin/env python

"""Tests of wrapper code in pycits."""

from pycits import blast
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal


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


def test_blastclust_execerr():
    """Correct error thrown if no blastclust"""
    try:
        bc = blast.Blastclust("blastclust_notexist")
    except NotExecutableError:
        return True
    else:
        return False
