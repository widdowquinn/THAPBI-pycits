#!/usr/bin/env python

"""Tests of wrapper code in pycits. vsearch"""

import os
import shutil
import subprocess

from pycits import vsearch
from pycits.tools import NotExecutableError, reformat_blast6_clusters

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data", "vsearch")
OUTDIR = os.path.join("tests", "test_out", "vsearch")
INFILE = os.path.join(INDIR, "test_db_for_abundance_checking.fasta")
DB = os.path.join(INDIR, "vsearch_tests_ITS_db.fasta")

# PARAMETERS
PREFIX = "test_run"
THRESHOLD = "0.96"
THREADS = "2"

# TARGET OUTPUT
TARGET_DEREP = os.path.join("tests", "test_targets", "vsearch",
                            "dedup_test_vsearch.fasta")
TARGET_BLAST6 = os.path.join("tests", "test_targets", "vsearch",
                             "target.blast6")
TARGET_UC_CLUSTERS = os.path.join("tests", "test_targets", "vsearch",
                                  "target.uc")
TARGET_R_FORMAT = os.path.join("tests", "test_targets", "vsearch",
                               "target.formatR")
TARGET_C_FAST_B6 = os.path.join("tests", "test_targets", "vsearch",
                                "clusterfast.blast6")
TARGET_C_FAST_UC = os.path.join("tests", "test_targets", "vsearch",
                                "target_runfast.clusters.uc")
TARGET_CENTROIDS = os.path.join("tests", "test_targets", "vsearch",
                                "target.centroids.fasta")
TARGET_ALIGNED = os.path.join("tests", "test_targets", "vsearch",
                              "target.alignedclusters.fasta")
TARGET_CONSENSUS = os.path.join("tests", "test_targets", "vsearch",
                                "target.consensus_cls_seq.fasta")


def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_vsearch_path():
    """vsearch is in $PATH"""
    vsearch = vsearch.Vsearch_derep("vsearch")





    
    
def get_sorted_list(in_file):
    """funct to return a sorted list.
    takes in a file, returns sorted list"""
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


def get_sorted_fa(fasta):
    """funct to return a sorted list of fa records.
    takes in a fasta, returns sorted list"""
    out_list = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        data = "%s\t%s" % (seq_record.id, seq_record.seq)
        out_list.append(data)
    return sorted(out_list)









@nottest    
def test_vsearch_cmd():
    """vsearch returns correct form of cmd-line"""
    outfname = os.path.join(OUTDIR, PREFIX +
                            'derep.fasta')
    derep = vsearch.Vsearch_derep("vsearch")
    target = ' '.join(["vsearch",
                       "--derep_fulllength", INFILE,
                       "--output", outfname,
                       "--sizeout"])
    assert_equal(derep.run(INFILE, OUTDIR, PREFIX,
                           dry_run=True), target)

@nottest    
def test_vsearch_exec_notexist():
    """error thrown if vsearch executable does not exist"""
    try:
        derep = vsearch.Vsearch_derep(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_vsearch_notexec():
    """Error thrown if vsearch not executable"""
    try:
        cluster = vsearch.Vsearch_derep("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_vsearch_exec():
    """Run vsearch on test data

    TODO: finer option could be passed??
    """
    derep = vsearch.Vsearch_derep("vsearch")
    # results are defined in the class:
    # factory class for vsearch class returned values
    # Results = namedtuple("Results", "command fastaout " +
    # "stdout stderr")
    # fasta_in, outdir, prefix, dry_run=False)
    result = derep.run(INFILE, OUTDIR, PREFIX)
    # use the named tuple to get the cluster results file
    vsearch_fasta = result.fasta

    with open(TARGET_DEREP, "rt") as target_fh:
        with open(vsearch_fasta, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


###################################################################
# Now to tests Vsearch_cluster
INFILE_CLUSTER = TARGET_DEREP

@nottest    
def test_Vsearch_cluster():
    """Vsearch_cluster instantiates with cmd-line if cd-hit is in $PATH"""
    cluster = vsearch.Vsearch_cluster("vsearch")

@nottest    
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

@nottest    
def test_Vsearch_cluster_exec_notexist():
    """Error thrown if Vsearch_cluster executable does not exist"""
    try:
        cluster = vsearch.Vsearch_cluster(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_Vsearch_cluster_notexec():
    """Error thrown if vsearch not executable"""
    try:
        cluster = vsearch.Vsearch_cluster("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_vsearch_exec():
    """Run vsearch on test data

    TODO: finer option could be passed??
    """
    cluster = vsearch.Vsearch_cluster("vsearch")
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


###################################################################
# Now to test conversion to another format
@nottest    
def test_convert_vsearch_format():
    """ testing function in tools to convert blast6 format to format
    for R"""
    cat_out = os.path.join(OUTDIR, "db_and_reads.fasta")
    cat_cmd = ' '.join(["cat",
                        DB,
                        INFILE,
                        ">",
                        cat_out])
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

    # reformat_blast6_clusters(blast6, db_and_reads, outfile)
    # call the function
    reformat_blast6_clusters(TARGET_BLAST6, cat_out,
                             os.path.join(OUTDIR, "tests.Rformat"))
    # convert to sorted lists
    result_R = get_sorted_list(os.path.join(OUTDIR, "tests.Rformat"))
    target_R = get_sorted_list(TARGET_R_FORMAT)
    assert_equal(result_R, target_R)


###################################################################
# Now to tests Vsearch_fastas

@nottest    
def test_vsearch_fastas():
    """Vsearch_fastas instantiates with cmd-line if cd-hit is in $PATH"""
    cluster = vsearch.Vsearch_fastas("vsearch")

@nottest    
def test_vsearch_fastas_cmd():
    """Vsearch_fastas instantiates and returns correct form of cmd-line"""
    cluster = vsearch.Vsearch_fastas("vsearch")
    target = ' '.join(["vsearch",
                       "--cluster_fast",
                       INFILE_CLUSTER,
                       "--id",
                       THRESHOLD,
                       "--centroids",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.centroids.fasta'),
                       "--msaout",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.alignedclusters.fasta'),
                       "--uc",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.fast.clusters.uc'),
                       "--consout",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.consensus_cls_seq.fasta'),
                       "--db",
                       DB,
                       "--threads",
                       THREADS,
                       "--blast6out",
                       os.path.join(OUTDIR,
                                    PREFIX +
                                    THRESHOLD +
                                    '.clusterfast.blast6')])
    # fasta_in, outdir, prefix, db, threads, threshold)
    assert_equal(cluster.run(INFILE_CLUSTER,
                             OUTDIR,
                             PREFIX,
                             DB,
                             THREADS,
                             THRESHOLD,
                             dry_run=True), target)

@nottest    
def test_vsearch_fastas_exec_notexist():
    """Error thrown if Vsearch_fastas executable does not exist"""
    try:
        cluster = vsearch.Vsearch_fastas(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_vsearch_fastas_notexec():
    """Error thrown if vsearch not executable"""
    try:
        cluster = vsearch.Vsearch_fastas("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False

@nottest    
def test_vsearch_exec():
    """Run vsearch on test data

    TODO: finer option could be passed??
    """
    cluster = vsearch.Vsearch_fastas("vsearch")

    # results are defined in the class:
    # factory class for vsearch class returned values
    # Results = namedtuple("Results", "command fastaout " +
    # "stdout stderr")
    # fasta_in, outdir, prefix, db, threads, threshold=0.99
    result2 = cluster.run(INFILE_CLUSTER,
                          OUTDIR,
                          PREFIX,
                          DB,
                          THREADS,
                          THRESHOLD)
    # use the named tuple to get the clustering results file
    # compare the blast6 output, convert them to sorted lists first
    target_blast = get_sorted_list(TARGET_C_FAST_B6)
    result_blast = get_sorted_list(result2.blast6)
    assert_equal(result_blast, target_blast)

    # UC outfiles
    target_uc = get_sorted_list(TARGET_C_FAST_UC)
    result_uc = get_sorted_list(result2.uc_clusters)
    assert_equal(result_uc, target_uc)

    # centroids outfiles
    target_cent = get_sorted_fa(TARGET_CENTROIDS)
    result_cent = get_sorted_fa(result2.centroids)
    assert_equal(target_cent, result_cent)

    # TARGET_ALIGNED
    target_alig = get_sorted_fa(TARGET_ALIGNED)
    result_alig = get_sorted_fa(result2.aligned)
    assert_equal(target_alig, result_alig)

    # TARGET_CONSENSUS
    target_conc = get_sorted_fa(TARGET_CONSENSUS)
    result_conc = get_sorted_fa(result2.consensus_cls)
    assert_equal(target_conc, result_conc)
