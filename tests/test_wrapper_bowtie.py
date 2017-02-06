#!/usr/bin/env python

"""Tests pycits wrapper for bowtie2.
Bowtie2 is supposed to be better for longer reads than bowtie1"""

import os
import shutil
from pycits import bowtie_build, bowtie_map
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
# bowtie2 wont work with .gz files. so need to decompress them.
import subprocess

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_bowtie")
DATABASE = os.path.join("data",
                        "ITS_db_NOT_conf_correct_last14bp_removd.fasta")
ASSEM_READS = os.path.join(INDIR, "pear_test.assembled.fastq.gz")
ASSEM_FA = os.path.join(INDIR, "bowtie",
                        "DNAMIX_S95_L001.assembled.fastq.bio.chopped.fasta.gz")

FA_INDEX = os.path.join("tests", "test_out_bowtie", "FA_INDEX")
THREADS = "1"

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "bowtie",
                      "DNAMIX_S95_L001_R1_001_bowtie.html")


def test_FQ_FA_is_compressed():
    """function to check the input data is already compressed"""
    file_list = [ASSEM_READS, ASSEM_FA]
    for i in file_list:
        if not os.path.isfile(i):
            cmd = comp(ASSEM_READS.split(".gz")[0])
            pipe = subprocess.run(decompcmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
        else:
            continue
    return True


def decomp(infname):
    """decompress a file"""
    cmd = ["gunzip",
           infname]
    cmd = ' '.join(cmd)
    return cmd


def comp(infname):
    """compress a file"""
    cmd = ["gzip",
           infname]
    cmd = ' '.join(cmd)
    return cmd


def test_outfolder():
    """checking to see if outfolder is made."""
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)


# test for indexing first
def test_bowtie():
    """bowtie_index instantiates with cmd-line if bowtie is in $PATH"""
    bowtie_build.Bowtie2_Build("bowtie2-build")


def test_bowtie_index_cmd():
    """bowtie2 instantiates, runs and returns correct form of cmd-line"""
    obj = bowtie_build.Bowtie2_Build("bowtie2-build")
    target = ' '.join(["bowtie2-build",
                       "--quiet",
                       "-f",
                       DATABASE,
                       FA_INDEX])
    assert_equal(obj.run(DATABASE, FA_INDEX, dry_run=True),
                 target)


def test_bowtie_exec_notexist():
    """Error thrown if bowtie_index executable does not exist"""
    try:
        obj = bowtie_build.Bowtie2_Build(os.path.join(".", "bowtie2-build"))
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie_notexec():
    """Error thrown if bowtie_index exe not executable"""
    try:
        obj = bowtie_build.Bowtie2_Build("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie_index_exec():
    """Run bowtie_index on test FA data and compare output
    to precomputed target"""
    obj = bowtie_build.Bowtie2_Build("bowtie2-build")
    result = obj.run(DATABASE, FA_INDEX)
    # test to see if it has produced the bai file.
    if not os.path.isfile(result.index + "1"):
        return False


def test_bowtie_index_exec():
    """Run bowtie_index on test FQ data and compare output
    to precomputed target"""
    obj = bowtie_build.Bowtie2_Build("bowtie2-build")
    result = obj.run(DATABASE, FA_INDEX)
    # test to see if it has produced the bai file.
    if not os.path.isfile(result.index + "1"):
        return False


#####################################################################
# now to test the mapping Bowtie2 scripts
# First using fasta files
def test_bowtie():
    """bowtie_map instantiates with cmd-line if bowtie is in $PATH"""
    bowtie_map.Bowtie2_Map("bowtie2")


def test_bowtie_map_cmd_using_fasta():
    """bowtie2 instantiates, runs and returns correct form of cmd-line
    fasta file input"""
    obj = bowtie_map.Bowtie2_Map("bowtie2")
    first = os.path.split(ASSEM_READS)[-1].split(".f")[0]
    second = "_Vs_" + os.path.split(FA_INDEX)[-1] + ".sam"
    outfile_name = first + second
    target = ' '.join(["bowtie2",
                       "--very-sensitive",
                       "--no-unal",
                       "-p", THREADS,
                       "-x",
                       FA_INDEX,
                       "-f",
                       ASSEM_FA.split(".gz")[0],
                       "-S",
                       os.path.join(OUTDIR, outfile_name)])
    assert_equal(obj.run(ASSEM_FA.split(".gz")[0], FA_INDEX,
                         OUTDIR, THREADS, dry_run=True),
                 target)


def test_bowtie_map_exec_notexist():
    """Error thrown if bowtie_map executable does not exist"""
    try:
        obj = bowtie_map.Bowtie2_Map(os.path.join(".", "bowtie2"))
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie_notexec():
    """Error thrown if bowtie_map exe not executable"""
    try:
        obj = bowtie_map.Bowtie2_Map("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie_map_cmd_fq():
    """bowtie2 map instantiates FQ input, runs and returns correct
    form of cmd-line"""
    obj = bowtie_map.Bowtie2_Map("bowtie2")
    first = os.path.split(ASSEM_READS)[-1].split(".f")[0]
    second = "_Vs_" + os.path.split(FA_INDEX)[-1] + ".sam"
    outfile_name = first + second
    target = ' '.join(["bowtie2",
                       "--very-sensitive",
                       "--no-unal",
                       "-p", THREADS,
                       "-x",
                       FA_INDEX,
                       "-q",
                       ASSEM_READS,
                       "-S",
                       os.path.join(OUTDIR, outfile_name)])
    assert_equal(obj.run(ASSEM_READS, FA_INDEX, OUTDIR, THREADS,
                 dry_run=True),
                 target)


def test_bowtie_map_FQ_exec():
    """Run bowtie_map on test FQ data and compare output
    to precomputed target"""
    # wont work with .gz files, so have to decomp these
    decompcmd = decomp(ASSEM_READS)
    pipe = subprocess.run(decompcmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    obj = bowtie_map.Bowtie2_Map("bowtie2")
    result = obj.run(ASSEM_READS.split(".gz")[0], FA_INDEX,
                     OUTDIR, THREADS)
    # compress them again
    comp_cmd = comp(ASSEM_READS.split(".gz")[0])
    pipe = subprocess.run(comp_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    # test to see if it has produced the sam file.
    if not os.path.isfile(result.sam):
        return False


def test_bowtie_map_FA_exec():
    """Run bowtie_map on test FA data and compare output
    to precomputed target"""
    decompcmd = decomp(ASSEM_FA)
    pipe = subprocess.run(decompcmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

    obj = bowtie_map.Bowtie2_Map("bowtie2")
    result = obj.run(ASSEM_FA.split(".gz")[0], FA_INDEX,
                     OUTDIR, THREADS)
    # compress them again
    comp_cmd = comp(ASSEM_FA.split(".gz")[0])
    pipe = subprocess.run(comp_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    # test to see if it has produced the sam file.

    if not os.path.isfile(result.sam):
        return False
