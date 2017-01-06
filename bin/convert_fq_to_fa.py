#!/usr/bin/env python
# coding: utf-8
# Title:
# script to reformat fastq to fasta format
# Why: after pear/ flash have assembled the reads.
# Need to convert from fastq to fa
# author: Peter Thorpe and Leighton Pritchard
# November 2016. The James Hutton Insitute, Dundee, UK.
# THIS HAS BEEN LEFT AS A STAND ALONE TOOL AS IT CAN BE USEFUL
# FOR OTHER APPLICATIONS


from Bio import SeqIO
import os
from sys import stdin, argv
import sys
import gzip
from optparse import OptionParser


def convert_fq_to_fa(in_file, out_file):
    """ Function to convert fq to fa
    -i the fastq file to be converted
    -o the desired out fasta file name
    the function will auto detect if the file is .gz
    requires biopython
    """
    # if the file is .gz we need to deal with it
    if in_file.endswith(".gz"):
        in_file = gzip.open(in_file, mode='rt', compresslevel=9,
                            encoding=None, errors=None,
                            newline=None)
    SeqIO.convert(in_file, "fastq", out_file, "fasta")

if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ convert_fq_to_fa -i in.fastq(.gz) -o outfile.fasta

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="fastq file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename, fasta",
                  metavar="FILE")


(options, args) = parser.parse_args()
in_file = options.in_file
out_file = options.out_file
(options, args) = parser.parse_args()

# run the program
convert_fq_to_fa(in_file, out_file)
