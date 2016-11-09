#!/usr/bin/env python

"""Tests of wrapper code in pycits. deduplicate"""

import os
import shutil

from pycits import deduplicate
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_deduplicate")
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)


def test_dedup():
    """deduplicate instantiates with cmd-line if deduplicate
    is in relative $PATH"""
    dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
    dedup = deduplicate.Deduplicate(dedup_prog)


def test_dedup_cmd():
    """deduplicate instantiates and returns correct form
of cmd-line"""
    dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
    test_fasta = os.path.join("tests", "test_data",
                             "dedup_test.fasta")
    outdirname = os.path.join("tests", "test_out_deduplicate")
    database = os.path.join(outdirname,
                            "test_database_out")
    dedup = deduplicate.Deduplicate(dedup_prog)
    out_file = os.path.join(outdirname,
                            "deduplicate_test.fasta")
    command =  ["python",
               dedup_prog,
               "--fasta", test_fasta,
               "--database", database,
               "-o", out_file]
    target = ' '.join(command)
    print (target)
    assert_equal(dedup.run(dedup_prog, test_fasta, database,
                           out_file), target)
    dedup.run(dedup_prog, test_fasta, database,
                           out_file)

test_dedup_cmd()

def test_dedup_exec_notexist():
    """Error thrown if dedup executable does not exist"""
    try:
        dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
        dedup = deduplicate.Deduplicate(os.path.join(".", dedup_prog))
    except NotExecutableError:
        return True
    else:
        return False


def test_dedup_notexec():
    """Error thrown if deduplicate dedup not executable"""
    try:
        dedup = deduplicate.Deduplicate("LICENSE")
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


def test_deduplicate_exec():
    """Run dedup on test data"""
    dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
    dedup = deduplicate.Deduplicate(dedup_prog)
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    
    outdirname = os.path.join("tests", "test_out_deduplicate")
    
    test_fasta = os.path.join("tests", "test_data",
                             "dedup_test.fasta")
    database = os.path.join("tests", "test_out_deduplicate",
                            "test_database_out")
    dedup = deduplicate.Deduplicate(dedup_prog)
    out_file = os.path.join("tests", "test_out_deduplicate",
                            "deduplicate_test.fasta")
    

    target_fasta = os.path.join("tests", "test_targets",
                          "dedup_test_coded.fasta")
    target_db = os.path.join("tests", "test_targets",
                          "test_database.out")
    result_fasta = os.path.join("tests", "test_out_deduplicate",
                                "deduplicate_test.fasta")
    result_db = os.path.join("tests", "test_out_deduplicate",
                             "test_database_out")

    dedup.run(dedup_prog, test_fasta, database,
                           out_file)    
    with open(target_fasta, "r") as target_fh:
        with open(result_fasta, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
            
    with open(target_db, "r") as target_fh:
        with open(result_db, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
