#!/usr/bin/env python3
#
# metapy_tools.py
#
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import gzip
from pycits.tools import convert_fq_to_fa
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# TODO: The PSL has the gzip library for this!!!!!
#       https://docs.python.org/3/library/gzip.html#examples-of-usage
def decompress(infile):
    """function to decompress gzipped reads"""
    cmd = ' '.join(["gunzip", infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    return infile.split(".gz")[0]  # os.splitext(infile)[0]!!!!!


def compress(infile):
    """function to compress reads, make them .gz"""
    cmd = ' '.join(["gzip", infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

# Report last exception as string


def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type,
                                              exc_value,
                                              exc_traceback))


def metapy_trim_seq(infname, outfname, lclip=53, rclip=0, minlen=100):
    """Trims all FASTA sequences in infname by lclip and rclip on the
    'left'- and 'right'-hand ends, respectively and writes the results to
    outfname. Any clipped sequence with length < minlen is not written.

    Defaults are equivalent to the default settings are for the current ITS1
    phy primers and the development ITS1 Phy database. PLEASE DONT BLINDLY
    USE THESE. The filenames, we now provide directly.
    """
    # ensure these are int, as they may have been passed as a string.
    lclip = int(lclip)
    rclip = int(rclip)
    with open(infname, 'r') as fh:
        # rclip needs to be set up to allow for zero values, using logic:
        # rclip_coordinate = len(s) - rclip
        # Use generators to save memory
        s_trim = (s[lclip:(len(s) - rclip)] for s in SeqIO.parse(fh, 'fasta'))
        return SeqIO.write((s for s in s_trim if len(s) >= minlen),
                           outfname, 'fasta')
