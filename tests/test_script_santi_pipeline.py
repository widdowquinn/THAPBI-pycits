#!/usr/bin/env python3

"""Tests of run_santi_pipeline.py script."""

import os
import shutil

from argparse import Namespace

from pycits.scripts import santi_pipeline

NAMESPACE = Namespace(blastclust='blastclust',
                      convert_format='convert_format',
                      fastqc='fastqc', indirname='data/reads',
                      join_paired_ends='join_paired_ends.py',
                      logfile='nosetest_santi_script_test.log',
                      muscle='muscle',
                      outdirname='nosetest_santi_script_test',
                      pick_closed_reference_otus='pick_closed_reference_otus.py',
                      pick_otus='pick_otus.py',
                      prefix='DNAMIX_S95_L001',
                      reference_fasta='data/database.fasta',
                      threads=4,
                      trim_quality='trim_quality',
                      verbose=True)


def test_complete_script():
    """complete script runs"""
    if os.path.isdir(NAMESPACE.outdirname):
        shutil.rmtree(NAMESPACE.outdirname)
    santi_pipeline.run_pycits_main(NAMESPACE)
