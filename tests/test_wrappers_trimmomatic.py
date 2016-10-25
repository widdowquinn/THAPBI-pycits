#!/usr/bin/env python

"""Tests of wrapper code in pycits. Trimmomatic"""

import os
import shutil

from pycits import trimmomatic
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

trimmo_prog = "~/scratch/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar"

INDIR = os.path.join("tests", "test_data")
outdir = os.path.join("tests", "test_out_trimmomatic")

# setting up some test variables
threads = 4
left_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq")
right_reads = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq")
read_prefix = left_reads.split("_R")[0]
read_prefix = read_prefix.split("/")[-1]
Left_outfile = os.path.join(outdir,
                    prefix+ "_paired_R1.fq.gz")
Right_outfile = os.path.join(outdir,
                    prefix+ "_paired_R2.fq.gz")
HEADCROP = 0

# This may not be suitbale for this program
#def test_trimmomatic():
#    """trimmomatic instantiates with cmd-line if trimmomatic is in $PATH"""
#    trim = trimmomatic.Trimmomatic("trimmomatic")


# test with 0 headcrop
def test_trimmomatic_cmd():
    """trimmomatic instantiates and returns correct form of cmd-line"""
    trim = trimmomatic.Trimmomatic(trimmo_path)
    HEADCROP = "HEADCROP:%s" %(str(HEADCROP))

    target = ' '.join(["java", "-jar",
               trimmo_prog,
               "PE",
               "-threads", str(threads),
               "-phred33",
               L_reads, R_reads,
               self._Left_outfile,
               "unpaired_R1.fq.gz",
               self._Right_outfile,
               "unpaired_R2.fq.gz",
               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3",
               HEADCROP,
               "TRAILING:3", "SLIDINGWINDOW:4:25", "MINLEN:70"])
    assert_equal(trim.run(trimmo_prog, left_reads, right_reads,
                 threads, "trimmed_reads", 0),
                 target)

# test with 40 headcrop
HEADCROP = "40"
def test_trimmomatic_cmd():
    """trimmomatic instantiates and returns correct form of cmd-line"""
    trim = trimmomatic.Trimmomatic(trimmo_prog)
    HEADCROP = "HEADCROP:%s" %(str(HEADCROP))
    Left_outfile = os.path.join(outdir,
                            prefix+ "_paired_R1.fq.gz")
    Right_outfile = os.path.join(outdir,
                            prefix+ "_paired_R2.fq.gz")
    target = ' '.join(["java", "-jar",
               trimmo_prog,
               "PE",
               "-threads", str(threads),
               "-phred33",
               L_reads, R_reads,
               Left_outfile,
               "unpaired_R1.fq.gz",
               Right_outfile,
               "unpaired_R2.fq.gz",
               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3",
               HEADCROP,
               "TRAILING:3", "SLIDINGWINDOW:4:25", "MINLEN:70"])
    assert_equal(trim.run(trimmo_prog, left_reads, right_reads,
                 threads, "trimmed_reads", 40),
                 target)

                 

### This may not be suitbale for this program
##def test_trimmomatic_exec_notexist():
##    """Error thrown if executable does not exist"""
##    try:
##        trim = blast.Blastclust(os.path.join(".",
##                    "~/scratch/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar"))
##    except NotExecutableError:
##        return True
##    else:
##        return False



def test_trimmomatic_exec():
    """Run trimmomatic on test data"""
    trim = trimmomatic.Trimmomatic(trimmo_prog)
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    trim.run((trimmo_prog, left_reads, right_reads,
                 threads, "trimmed_reads", 40)

    left_result = Left_outfile
    right_result = Right_outfile
    left_test = os.path.join("tests", "test_targets",
                          "DNAMIX_S95_L001_paired_R1.fq.gz")
    right_test = os.path.join("tests", "test_targets",
                          "DNAMIX_S95_L001_paired_R2.fq.gz")
    with open(left_test, "r") as left_test_fh:
        with open(left_result, "r") as left_result_fh:
            assert_equal(left_result_fh.read(), left_test_fh()
                         
    with open(right_test, "r") as right_test_fh:
        with open(right_result, "r") as right_result_fh:
            assert_equal(right_result_fh.read(), right_test_fh.read())
                      
