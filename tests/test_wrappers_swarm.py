#!/usr/bin/env python

"""Tests of wrapper code in pycits. swarm"""

import os
import shutil

from pycits import swarm
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_swarm")
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)


def test_dedup():
    """swarm instantiates with cmd-line if swarm
    is in $PATH"""
    cluster = swarm.Swarm("swarm")


def test_dedup_cmd():
    """swarm instantiates and returns correct form
    of cmd-line"""
    test_fasta = os.path.join("tests", "test_data",
                             "swarm_coded_with_abundance.fasta")
    outdirname = os.path.join("tests", "test_out_swarm")
    #set up the class
    cluster = swarm.Swarm("swarm")
    clustering_threshold = 1
    clustering_threshold = str(clustering_threshold)
    out_file = os.path.join(outdirname,
                            "swarm_clustering_d%s" %
                            (clustering_threshold),
                            "swarm_clustering_d%s" %
                            (clustering_threshold))
    out = os.path.join(outdirname,
                       "swarm_clustering_d%s" %
                       (clustering_threshold))
    command =  ["swarm",
                "-t",
                "1",
                "-d",
                "1",
                "-o",
                out_file,
                test_fasta]
    target = ' '.join(command)
    print ("\n\ntarget = ",target)
    class_dir, class_command = cluster.run(test_fasta, "1",\
                                           1, outdirname)
    print ("\n\nclass_command = ", class_command)
    assert_equal(class_command, target)
    cluster.run(test_fasta, "1", 1, outdirname)

test_dedup_cmd()

def test_dedup_exec_notexist():
    """Error thrown if dedup executable does not exist"""
    try:
        cluster = swarm.Swarm(os.path.join(".", "swarm"))
    except NotExecutableError:
        return True
    else:
        return False


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
    #set up the class
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

