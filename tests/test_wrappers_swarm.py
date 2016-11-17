#!/usr/bin/env python

"""Tests of wrapper code in pycits. swarm"""

import os
import shutil

from pycits import swarm
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "swarm")
OUTDIR = os.path.join("tests", "test_out", "swarm")
INFILE = os.path.join(INDIR, "swarm_coded_with_abundance.fasta")
OUTFILE = os.path.join(OUTDIR, "swarm.out")

def test_swarm():
    """swarm instantiates with cmd-line if swarm is in $PATH"""
    cluster = swarm.Swarm("swarm")


def test_swarm_cmd():
    """swarm instantiates and returns correct form of cmd-line"""
    cluster = swarm.Swarm("swarm")
    target = ' '.join(["swarm -t 1 -d 1",
                       "-o {0}".format(OUTFILE),
                       INFILE])
    assert_equal(cluster.run(INFILE, OUTDIR, 1, 1, dry_run=True), target)


@nottest
def test_dedup_exec_notexist():
    """Error thrown if dedup executable does not exist"""
    try:
        cluster = swarm.Swarm(os.path.join(".", "swarm"))
    except NotExecutableError:
        return True
    else:
        return False


@nottest
def test_dedup_notexec():
    """Error thrown if swarm dedup not executable"""
    try:
        cluster = swarm.Swarm("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def build_diff_cmd(infname1, infname2):
    """Build a command-line for diff"""
    cmd = ["diff",
           infname1,
           infname2]
    build_diff_cmd = ' '.join(cmd)
    return build_diff_cmd


@nottest
def test_swarm_exec():
    """Run dedup on test data"""
    cluster = swarm.Swarm("swarm")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

    test_fasta = os.path.join("tests", "test_data",
                              "swarm_coded_with_abundance.fasta")
    outdirname = os.path.join("tests", "test_out_swarm")
    # set up the class
    clustering_threshold = 1
    clustering_threshold = str(clustering_threshold)
    out_file = os.path.join(outdirname,
                            "swarm_clustering_d%s" %
                            (clustering_threshold))

    target = os.path.join("tests", "test_targets",
                          "swarm_test_d1.out")
    result = os.path.join("tests", "test_out_swarm",
                          "swarm_clustering_d1", "swarm_clustering_d1")
    cluster.run(test_fasta, "1", 1, outdirname)
    with open(target, "r") as target_fh:
        with open(result, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
