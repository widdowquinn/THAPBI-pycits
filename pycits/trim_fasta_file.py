#!/usr/bin/env python
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe
# os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

# Title:
# script to trim sequences to the region of ITS/
# or the given region, as defined by the user.

########################################################
# functions

def check_db_abundance(title):
    """function to check that the database entries have
    abundance _1 appended to the end of the title for
    swarm. If they dont it will add them"""
    if not title.endswith("_1"):
        title += "_1"
    return title
    

def reformat_fasta_name(fasta, database, 
                        left_trim, 
                        right_trim, out):
    """function to retun the fa file catatentated with
    the databse seq. The assembled reads will be trimmed
    at the given lengths. """
    f= open(out, 'w')
    left_trim = int(left_trim)
    right_trim = int(right_trim)
    # iterate through the database. append abundance
    for seq_record in SeqIO.parse(database, "fasta"):
        seq_record.id = check_db_abundance(seq_record.id)
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")
    # iterate through the assembled file. Trim the seq
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if left_trim or right_trim:
            #if given left and right primer. Slice these off
            seq_record.seq = seq_record.seq[left_trim:\
                                            len(seq_record.seq) - right_trim]
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")
    f.close()



###########################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python trim_fasta_file.py -f seq.fasta --database ITS.fasta
    --left 53 --right 0

--left is the length of the seq to chop off, from left
        - 0 by default
--right as above but from the right = 0 by default


requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered",
                  metavar="FILE")

parser.add_option("-d","--database", dest="database",
                  default=None,
                  help="fasta file to have names altered",
                  metavar="FILE")

parser.add_option("-r", "--right", dest="right_trim", default=0,
                  help="length right primer "
                  "or the amount of seq to trim from the right")

parser.add_option("-l", "--left", dest="left_trim", default=0,
                  help="length of left primer "
                  "or the amount to trim from the left. If using "
                  "Phytophthora ITS1 primers set this to 53")

parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filenames")


(options, args) = parser.parse_args()

fasta = options.fasta
database = options.database
left_trim = options.left_trim
right_trim = options.right_trim
out = options.out

# run the program

# print "ITS = ", ITS

# biopython imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_name(fasta, database,
                    left_trim, right_trim,
                    out)

