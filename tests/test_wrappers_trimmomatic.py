#!/usr/bin/env python

"""Tests of wrapper code in pycits. Trimmomatic"""

import os
import shutil
import sys

import subprocess

from pycits import trimmomatic
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal


def test_trimmomatic_cmd():
    """trimmomatic instantiates and returns correct form of cmd-line"""
    trimmo_path = "/home/pt40963/scratch/Downloads/" + \
        "Trimmomatic-0.36/trimmomatic-0.36.jar"

    INDIR = os.path.join("tests", "test_data")
    OUTDIR = os.path.join("tests", "test_out_trimmomatic")

    threads = 4
    left_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
    right_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")

    read_prefix = left_reads.split("_R")[0]
    read_prefix = read_prefix.split("/")[-1]

    Left_outfile = os.path.join(OUTDIR,
                                read_prefix + "_paired_R1.fq.gz")
    Right_outfile = os.path.join(OUTDIR,
                                 read_prefix + "_paired_R2.fq.gz")

    # test with 40 headcrop
    HEADCROP = "40"
    trim = trimmomatic.Trimmomatic(trimmo_path)
    HEADCROP = "HEADCROP:%s" % (str(HEADCROP))
    target = ' '.join(["java", "-jar",
                       trimmo_path,
                       "PE",
                       "-threads", str(threads),
                       "-phred33",
                       left_reads,
                       right_reads,
                       Left_outfile,
                       "unpaired_R1.fq.gz",
                       Right_outfile,
                       "unpaired_R2.fq.gz",
                       "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3",
                       HEADCROP,
                       "TRAILING:3", "SLIDINGWINDOW:4:25", "MINLEN:70"])
    trim.run(trimmo_path, left_reads, right_reads,
             threads, OUTDIR, 40)

# run the tests with headcrop 40
test_trimmomatic_cmd()


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


def test_trimmomatic_exec():
    """Run trimmomatic on test data"""
    trimmo_path = "/home/pt40963/scratch/Downloads/" + \
        "Trimmomatic-0.36/trimmomatic-0.36.jar"
    INDIR = os.path.join("tests", "test_data")
    OUTDIR = os.path.join("tests", "test_out_trimmomatic")
    threads = 4
    left_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
    right_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
    read_prefix = left_reads.split("_R")[0]
    read_prefix = read_prefix.split("/")[-1]
    Left_outfile = os.path.join(OUTDIR,
                                read_prefix + "_paired_R1.fq.gz")
    Right_outfile = os.path.join(OUTDIR,
                                 read_prefix + "_paired_R2.fq.gz")
    # test with 40 headcrop
    HEADCROP = "40"
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    trim = trimmomatic.Trimmomatic(trimmo_path)
    trim.run(trimmo_path, left_reads, right_reads,
             threads, "tests/test_out_trimmomatic", 40)

    left_test = os.path.join("tests", "test_targets",
                             "DNAMIX_S95_L001_paired_R1.md5sum")
    right_test = os.path.join("tests", "test_targets",
                              "DNAMIX_S95_L001_paired_R2.md5sum")

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

    left_diff_test = build_diff_cmd(left_test,
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
    right_diff_test = build_diff_cmd(right_test,
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


test_trimmomatic_exec()
