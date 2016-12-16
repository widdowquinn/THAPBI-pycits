#!/usr/bin/env python

"""Tests of wrapper code in pycits. spades.py for
error correction
"""

import os
import shutil
import subprocess
from pycits import error_correction
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_EC")
READS1_gz = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2_gz = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = "DNAMIX_S95_L001_"
# TARGET OUTPUT DATA
TARGET_L_gz = os.path.join("tests", "test_targets", "error_correction",
                      "DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq.gz")
TARGET_R_gz = os.path.join("tests", "test_targets", "error_correction",
                      "DNAMIX_S95_L001_R2_001.fastq.00.0_0.cor.fastq.gz")

def build_decompress(infname):
    """Build a command-line for md5sum"""
    cmd = ["gunzip",
           infname]
    cmd = ' '.join(cmd)
    return cmd

def build_compress(infname):
    """Build a command-line for md5sum"""
    cmd = ["gzip",
           infname]
    cmd = ' '.join(cmd)
    return cmd

# in order to compare the file for their content
# we need these deompressed.
# we will get Linux to do this for us.
read_to_decomp = [READS1_gz, READS2_gz, TARGET_L_gz, TARGET_R_gz]
for i in read_to_decomp:
    # check if the file exists, so idont have to keep compressing during dev
    if not os.path.isfile(i):
        continue
    print ("decomp :", i)
    decomp_command = build_decompress(i)
    subprocess.run(decomp_command, shell=True,
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq")
TARGET_L = os.path.join("tests", "test_targets", "error_correction",
                      "DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq")
TARGET_R = os.path.join("tests", "test_targets", "error_correction",
                      "DNAMIX_S95_L001_R2_001.fastq.00.0_0.cor.fastq")


def sort_fq_output(fq_file):
    """function to sort the output of the EC assembled data
    . Normally this is returned in an unordered manner. Needs
    sorting for tests.
    returned as a sorted(set)"""
    # open the file
    try:
        in_file = open(fq_file)
    except:
        ValueError
        fq_file = fq_file.split(".gz")[0]
        in_file = open(fq_file)
    fq_set = set([])
    for (title, seq, qual) in (FastqGeneralIterator(in_file)):
        read_info = "%s\n+\n%s\n" % (seq, qual)
        fq_set.add(title + read_info)
    in_file.close()
    return sorted(fq_set)
##########################################################################
def test_Error_Correction():
    """spades.py instantiates with cmd-line if spades.py is in $PATH"""
    ec = error_correction.Error_Correction("spades.py")


def test_Error_Correction_cmd():
    """spades.py instantiates, runs and returns correct form of cmd-line"""
    obj = error_correction.Error_Correction("spades.py")
    target = ' '.join(["spades.py",
                       "-1",
                       READS1,
                       "-2",
                       READS2,
                       "--only-error-correction",
                       "--threads", "4",
                       "-o", OUTDIR])
    assert_equal(obj.run(READS1, READS2, 4,
                         os.path.join(OUTDIR),
                         dry_run=True),target)


#######################################################################

def test_spade_py_exec_notexist():
    """Error thrown if spade.py executable does not exist"""
    try:
        ec = error_correction.Error_Correction(os.path.join(".", "spades.py"))
    except NotExecutableError:
        return True
    else:
        return False

def test_EC_notexec():
    """Error thrown if EC exe not executable"""
    try:
        obj = error_correction.Error_Correction("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_error_correction_exec():
    """Run spades.py on test data and compare output to precomputed target
    """
    obj = error_correction.Error_Correction("spades.py")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = obj.run(READS1, READS2, 4, OUTDIR, dry_run=False)
    #print ("results:", result.right_read_correct)
    # call function to sort the data. returned as a sorted(set)
    # the output is not sorted, so a direct comparason of the file
    # does not work. So have to sort first.
    #print ("now sorting TARGET_L", TARGET_L)
    targed_data_sorted_left = sort_fq_output(TARGET_L)
    #print ("now sorting TARGET_R")
    targed_data_sorted_right = sort_fq_output(TARGET_R)
    # Left_read_correct - name from named tuple and originally
    # from build_command
    # first need to decompress these. 
    read_to_decomp = [result.Left_read_correct, result.right_read_correct]
    #print("decompressing: ", read_to_decomp)
    for i in read_to_decomp:
        # check if the file exists, so idont have to keep compressing during dev
        if not os.path.isfile(i):
            continue
        decomp_command = build_decompress(i)
        subprocess.run(decomp_command, shell=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    test_data_sorted_L = sort_fq_output(result.Left_read_correct)
    test_data_sorted_R = sort_fq_output(result.right_read_correct)
    # test if they are equal
    assert_equal(targed_data_sorted_left, test_data_sorted_L)
    assert_equal(targed_data_sorted_right, test_data_sorted_R)


############################################################################
# To clean up compress the origin data file
# compress them  
for i in read_to_decomp:
    name = i.split(".gz")[0]
    comp_command = build_compress(name)
    subprocess.run(comp_command, shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True)




