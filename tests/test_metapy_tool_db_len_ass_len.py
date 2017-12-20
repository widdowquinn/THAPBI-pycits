#!/usr/bin/env python3

"""Test program for a function in the module metapy_tools:
db_len_assembled_len_ok
"""

import os
import shutil

from nose.tools import nottest, assert_almost_equal
from pycits.metapy_tools import db_len_assembled_len_ok
from pycits.tools import NotExecutableError

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "metapy_tools")
OUTDIR = os.path.join("tests", "test_out", "metapy_tools")
DB = os.path.join(INDIR, "testdb.fasta")
INFILE_TOO_SHORT = os.path.join(INDIR,
                                "assembled_fa_too_short.fasta")
INFILE_TOO_LONG = os.path.join(INDIR,
                               "assembled_fa_too_long.fasta")
INFILE_GOOD = os.path.join(INDIR,
                           "assembled_fa_good.fasta")

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "metapy_tools",
                      "passed_target.fasta")


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def seqfile_to_length_list(filename, fmt='fasta'):
    """Returns list of sequence lengths for passed filename."""
    return [len(seq) for seq in SeqIO.parse(filename, fmt)]


def test_assem_db_seq_len_short():
    """Run db_len_assem_len_ok too SHORT from metapy_tools
    """
    error, vals = db_len_assembled_len_ok(seqfile_to_length_list(DB),
                                          seqfile_to_length_list(INFILE_TOO_SHORT))
    (db_mean, db_sd, assemb_mean, assemb_sd) = vals
    assert error == "fail\t-"
    assert_almost_equal(db_mean, 196.953, delta=0.001)
    assert_almost_equal(db_sd, 18.846, delta=0.001)
    assert_almost_equal(assemb_mean, 21.636, delta=0.001)
    assert_almost_equal(assemb_sd, 8.391, delta=0.001)


def test_assem_db_seq_len_long():
    """Run db_len_assem_len_ok too LONG from metapy_tools
    """
    error, vals = db_len_assembled_len_ok(seqfile_to_length_list(DB),
                                          seqfile_to_length_list(INFILE_TOO_LONG))
    (db_mean, db_sd, assemb_mean, assemb_sd) = vals
    assert error == "fail\t+"
    assert_almost_equal(db_mean, 196.953, delta=0.001)
    assert_almost_equal(db_sd, 18.846, delta=0.001)
    assert_almost_equal(assemb_mean, 1150.8, delta=0.001)
    assert_almost_equal(assemb_sd, 141.426, delta=0.001)


def test_assem_db_seq_len_ok():
    """Run db_len_assem_len_ok OK from metapy_tools
    """
    error, vals = db_len_assembled_len_ok(seqfile_to_length_list(DB),
                                          seqfile_to_length_list(INFILE_GOOD),
                                          10)
    (db_mean, db_sd, assemb_mean, assemb_sd) = vals
    assert error == "ok"
    assert_almost_equal(db_mean, 196.953, delta=0.001)
    assert_almost_equal(db_sd, 18.846, delta=0.001)
    assert_almost_equal(assemb_mean, 182.684, delta=0.001)
    assert_almost_equal(assemb_sd, 19.262, delta=0.001)


def test_assem_db_seq_len_ok_sd_0():
    """Run db_len_assem_len 0 sd threshold
    """
    error, msg = db_len_assembled_len_ok(seqfile_to_length_list(DB),
                                         seqfile_to_length_list(INFILE_GOOD),
                                         0)
    assert error.startswith("fail")
