#!/usr/bin/env python

"""Tests of wrapper code in pycits. rename"""

import os
import shutil

from pycits import rename_clusters_new_to_old
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_rename")
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)


def test_rename():
    """rename instantiates with cmd-line if rename
    is in relative $PATH"""
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    rename = rename_clusters_new_to_old.Rename(rename_prog)


def test_rename_cmd():
    """rename instantiates and returns correct form
    of cmd-line"""
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    test_in = os.path.join("tests", "test_data",
                             "clustered_tests_for_renaming.out")
    outdirname = os.path.join("tests", "test_out_rename")
    database = os.path.join("tests", "test_data",
                            "database_for_renaming.out")
    rename = rename_clusters_new_to_old.Rename(rename_prog)
    out_file = os.path.join(outdirname,
                            "rename_test.txt")
    command =  ["python",
               rename_prog,
               "-i", test_in,
               "--database", database,
               "-o", out_file]
    target = ' '.join(command)
    print (target)
    assert_equal(rename.run(rename_prog, test_in, database,
                           out_file), target)
    rename.run(rename_prog, test_in, database,
                           out_file)

test_rename_cmd()

def test_rename_exec_notexist():
    """Error thrown if rename executable does not exist"""
    try:
        rename_prog = os.path.join("pycits",
                                   "parse_clusters_new_to_old_name.py")
        rename = rename_clusters_new_to_old.Rename(os.path.join
                                                  (".", rename_prog))
    except NotExecutableError:
        return True
    else:
        return False


def test_rename_notexec():
    """Error thrown if rename rename not executable"""
    try:
        rename = rename_clusters_new_to_old.Rename("LICENSE")
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


def test_rename_exec():
    """Run rename on test data"""
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    rename = rename_clusters_new_to_old.Rename(rename_prog)
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    test_in = os.path.join("tests", "test_data",
                           "clustered_tests_for_renaming.out")
    outdirname = os.path.join("tests", "test_out_rename")
    database = os.path.join("tests", "test_data",
                            "database_for_renaming.out")
    rename = rename_clusters_new_to_old.Rename(rename_prog)
    out_file = os.path.join(outdirname,
                            "rename_test.txt")
    

    target = os.path.join("tests", "test_targets",
                          "rename.result")

    result = out_file

    rename.run(rename_prog, test_in, database,
               out_file)  
    with open(target, "r") as target_fh:
        with open(result, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
            
