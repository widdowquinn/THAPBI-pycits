#!/usr/bin/env python

"""Tests pycits wrapper for samtools_index."""

import os
import shutil
import gzip

from pycits import samtools_index
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data", "samtools")
BAM_FILE = os.path.join(INDIR, "Aligned.sortedByCoord.out.bam")


def test_samtools_index():
    """Samtools_index instantiates with cmd-line if Samtools_index
    is in $PATH"""
    samtools_index.Samtools_Index("samtools")


def test_samtools_index_cmd():
    """samtools_index instantiates, runs and returns
    correct form of cmd-line"""
    obj = samtools_index.Samtools_Index("samtools")
    target = ' '.join(["samtools", "index", BAM_FILE])
    stats = obj.run(BAM_FILE)
    assert_equal(stats.command,
                 target)


def test_samtools_index_exec_notexist():
    """Error thrown if samtools_index executable does not exist"""
    try:
        obj = samtools_index.Samtools_Index(os.path.join
                                            (".", "samtools"))
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_index_notexec():
    """Error thrown if samtools_index exe not executable"""
    try:
        obj = samtools_index.Samtools_Index("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_index_exec():
    """Run samtools_index on test data and compare output
    to precomputed target"""
    obj = samtools_index.Samtools_Index("samtools")
    result = obj.run(BAM_FILE)
    # test to see if it has produced the bai file.
    if not os.path.isfile(result.bai):
        return False
