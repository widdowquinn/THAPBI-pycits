#!/usr/bin/env python

"""Tests of wrapper code in pycits. swarm"""

import os
import shutil

from pycits import swarm
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "swarm")
OUTDIR = os.path.join("tests", "test_out_swarm")
INFILE = os.path.join(INDIR, "swarm_coded_with_abundance.fasta")
OUTFILE = os.path.join(OUTDIR, "swarm.out")

# TARGET OUTPUT
TARGET = os.path.join("tests", "test_targets", "swarm", "swarm.out")

def test_swarm():
    """swarm instantiates with cmd-line if swarm is in $PATH"""
    cluster = swarm.Swarm("swarm")


def test_swarm_cmd():
    """swarm instantiates and returns correct form of cmd-line"""
    cluster = swarm.Swarm("swarm")
    target = ' '.join(["swarm -t 1 -d 1",
                       "-o {0}".format(OUTFILE),
                       INFILE])
    print(target)
    assert_equal(cluster.run(INFILE, OUTDIR, 1, 1, dry_run=True), target)


def test_swarm_exec_notexist():
    """Error thrown if swarm executable does not exist"""
    try:
        cluster = swarm.Swarm(os.path.join(".", "swarm"))
    except NotExecutableError:
        return True
    else:
        return False


def test_swarm_notexec():
    """Error thrown if swarm not executable"""
    try:
        cluster = swarm.Swarm("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_swarm_exec():
    """Run swarm on test data"""
    cluster = swarm.Swarm("swarm")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = cluster.run(INFILE, OUTDIR, 1, 1)
    with open(TARGET, "r") as target_fh:
        with open(result[0], "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
