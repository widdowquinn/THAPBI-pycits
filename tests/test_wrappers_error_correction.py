#!/usr/bin/env python

"""Tests of wrapper code in pycits. Trimmomatic"""

import os
import shutil
import subprocess

from pycits import error_correction
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_EC")

# this will depend on the binarie name of EC
# WE will need to change this!!


def test_EC():
    """EC instantiates with cmd-line if EC is in $PATH
    SPAdes"""
    ByHam = "/home/pt40963/scratch/Downloads/SPAdes-3.9.0-Linux/bin/spades.py"
    error_correct = error_correction.Error_correction(ByHam)


def test_EC_cmd():
    """EC instantiates and returns correct form of cmd-line"""
    ByHam = "/home/pt40963/scratch/Downloads/SPAdes-3.9.0-Linux/bin/spades.py"
    error_correct = error_correction.Error_correction(ByHam)
    command = ["python",
            ByHam,
            "-1",
            "tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
            "-2",
            "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
            "--only-error-correction",
            "--threads",
            "4",
            "-o",
            "tests/test_out_EC"]
    target = ' '.join(command)
    assert_equal(error_correct.run(ByHam ,
                                   "tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
                                   "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
                                   "tests/test_out_EC", "4"), target)


def test_EC_exec_notexist():
    """Error thrown if EC executable does not exist"""
    ByHam = "/home/pt40963/scratch/Downloads/SPAdes-3.9.0-Linux/bin/spades.py"
    try:
        error_correct = error_correction.Error_correction(os.path.join(".", ByHam))
    except NotExecutableError:
        return True
    else:
        return False


def test_EC_notexec():
    """Error thrown if EC exe not executable"""
    try:
        error_correct = error_correction.Error_correction("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def build_cmd(infname):
    """Build a command-line for md5sum"""
    cmd = ["md5sum",
           infname,
           "|",
           "cut",
           "-d",
           "' '",
           "-f1",
           ">",
           infname + ".md5sum"]
    cmd = ' '.join(cmd)
    return cmd


def build_diff_cmd(infname1, infname2):
    """Build a command-line for diff"""
    cmd = ["diff",
           infname1,
           infname2]
    build_diff_cmd = ' '.join(cmd)
    return build_diff_cmd


def test_EC_exec():
    """Run EC on test data"""
    ByHam = "/home/pt40963/scratch/Downloads/SPAdes-3.9.0-Linux/bin/spades.py"
    error_correct = error_correction.Error_correction(ByHam)
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    error_correct.run(ByHam,
                      "tests/test_data/DNAMIX_S95_L001_R1_001.fastq.gz",
                      "tests/test_data/DNAMIX_S95_L001_R2_001.fastq.gz",
                      "tests/test_out_EC", "4")
    L_target = os.path.join("tests", "test_targets",
                          "DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq.gz.md5sum")
    R_target = os.path.join("tests", "test_targets",
                          "DNAMIX_S95_L001_R2_001.fastq.00.0_0.cor.fastq.gz.md5sum")
    Left_outfile = os.path.join("tests", "test_out_EC", "corrected",
                          "DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq.gz")
    Right_outfile = os.path.join("tests", "test_out_EC", "corrected",
                          "DNAMIX_S95_L001_R2_001.fastq.00.0_0.cor.fastq.gz")

    # builds: md5sum Left_outfile | cut -d ' ' -f1 > Left_outfile.md5sum
    left_test_command = build_cmd(Left_outfile)
    right_test_command = build_cmd(Right_outfile)
    L_pipe = subprocess.run(left_test_command, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True)
    print ("\n\nMD5SUM : = ", L_pipe, "\n\n")
    if L_pipe.returncode != 0:
        print("md5sum fail on left file")
        sys.exit(1)
    R_pipe = subprocess.run(right_test_command, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True)
    if R_pipe.returncode != 0:
        print("md5sum fail on right file")
        sys.exit(1)
        
    # linx diff command to test the md5sums
    left_diff_test = build_diff_cmd(L_target,
                                    Left_outfile + ".md5sum")

    print ("left_diff_test = ", left_diff_test, "\n\n")
    L_D_pipe = subprocess.run(left_diff_test, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    print (L_D_pipe)
    if L_D_pipe.returncode != 0:
        print("left files do not match")
        sys.exit(1)
    right_diff_test = build_diff_cmd(R_target,
                                     Right_outfile + ".md5sum")
    R_D_pipe = subprocess.run(right_diff_test, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    pipe = subprocess.run(right_test_command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    print ("\n R_D_pipe : ", R_D_pipe)
    if L_D_pipe.returncode != 0:
        print("right files do not match")
        sys.exit(1)
