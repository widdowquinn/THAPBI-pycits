#!/usr/bin/env python

"""Tests of wrapper code in pycits. vsearch"""

import os
import shutil

from pycits import vsearch
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "vsearch")
OUTDIR = os.path.join("tests", "test_out_vsearch")
INFILE_DEREP = os.path.join(INDIR, "test_db_for_abundance_checking.fasta")
PREFIX = "test_run"
THRESHOLD = "0.96"
THREADS = "2"
DB = os.path.join("data", "ITS_db_NOT_conf_correct_last14bp_removd.fasta")

# TARGET OUTPUT
TARGET_DEREP = os.path.join("tests", "test_targets", "vsearch",
                            "dedup_test_vsearch.fasta")
TARGET_BLAST6 = os.path.join("tests", "test_targets", "vsearch",
                             "target.blast6")
TARGET_UC_CLUSTERS = os.path.join("tests", "test_targets", "vsearch",
                                  "target.uc")


def get_sorted_list(in_file):
    """funct to return a sorted list.
    takes in a file, restunr sorted list"""
    out_list = []
    with open(in_file) as fh:
        data = fh.read().split("\n")
        for line in data:
            if not line.strip():
                continue  # if the last line is blank
            if line.startswith("#"):  # dont want comment lines
                continue
        out_list.append(line)
    return sorted(out_list)


def test_Vsearch_derep():
    """Vsearch_derep instantiates with cmd-line if cd-hit is in $PATH"""
    derep = vsearch.Vsearch_derep("vsearch")


def test_Vsearch_derep_cmd():
    """Vsearch_derep instantiates and returns correct form of cmd-line"""
    outfname = os.path.join(OUTDIR, PREFIX +
                            'derep.fasta')
    derep = vsearch.Vsearch_derep("vsearch")
    target = ' '.join(["vsearch",
                       "--derep_fulllength", INFILE_DEREP,
                       "--output", outfname,
                       "--sizeout"])
    assert_equal(derep.run(INFILE_DEREP, OUTDIR, PREFIX,
                           dry_run=True), target)


def test_Vsearch_derep_exec_notexist():
    """Error thrown if Vsearch_derep executable does not exist"""
    try:
        derep = vsearch.Vsearch_derep(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False


def test_Vsearch_derep_notexec():
    """Error thrown if vsearch not executable"""
    try:
        cluster = vsearch.Vsearch_derep("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_vsearch_exec():
    """Run vsearch on test data

    TODO: finer option could be passed??
    """
    derep = vsearch.Vsearch_derep("vsearch")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

    # results are defined in the class:
    # factory class for vsearch class returned values
    # Results = namedtuple("Results", "command fastaout " +
    # "stdout stderr")
    # fasta_in, outdir, prefix, dry_run=False)
    result = derep.run(INFILE_DEREP, OUTDIR, PREFIX)
    # use the named tuple to get the cluster results file
    vsearch_fasta = result.fasta

    with open(TARGET_DEREP, "rt") as target_fh:
        with open(vsearch_fasta, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


###################################################################
# Now to tests Vsearch_cluster
INFILE_CLUSTER = TARGET_DEREP


def test_Vsearch_cluster():
    """Vsearch_cluster instantiates with cmd-line if cd-hit is in $PATH"""
    cluster = vsearch.Vsearch_cluster("vsearch")


def test_Vsearch_cluster_cmd():
    """Vsearch_cluster instantiates and returns correct form of cmd-line"""
    cluster = vsearch.Vsearch_cluster("vsearch")
    target = ' '.join(["vsearch",
                       "--usearch_global",
                       INFILE_CLUSTER,
                       "--id",
                       THRESHOLD,
                       "--uc",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.clusters.uc'),
                       "--db",
                       DB,
                       "--threads",
                       str(THREADS),
                       "--blast6out",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD + '.blast6')])
    assert_equal(cluster.run(INFILE_CLUSTER,
                             OUTDIR,
                             PREFIX,
                             DB,
                             THREADS,
                             THRESHOLD,
                             dry_run=True), target)


def test_Vsearch_cluster_exec_notexist():
    """Error thrown if Vsearch_cluster executable does not exist"""
    try:
        cluster = vsearch.Vsearch_cluster(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False


def test_Vsearch_cluster_notexec():
    """Error thrown if vsearch not executable"""
    try:
        cluster = vsearch.Vsearch_cluster("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_vsearch_exec():
    """Run vsearch on test data

    TODO: finer option could be passed??
    """
    cluster = vsearch.Vsearch_cluster("vsearch")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

    # results are defined in the class:
    # factory class for vsearch class returned values
    # Results = namedtuple("Results", "command fastaout " +
    # "stdout stderr")
    # fasta_in, outdir, prefix, db, threads, threshold=0.99
    result = cluster.run(INFILE_CLUSTER,
                         OUTDIR,
                         PREFIX,
                         DB,
                         THREADS,
                         THRESHOLD)
    # use the named tuple to get the clustering results file
    # compare the blast6 output, convert them to sorted lists first
    target_blast = get_sorted_list(TARGET_BLAST6)
    result_blast = get_sorted_list(result.blast6)
    assert_equal(result_blast, target_blast)

    # UC outfiles
    target_uc = get_sorted_list(TARGET_UC_CLUSTERS)
    result_us = get_sorted_list(result.uc_clusters)
    assert_equal(result_us, target_uc)
