#!/usr/bin/env python
# Title:
# script to reformt fasta ID names and append the abundance of
# that seq at the end of the new name.
# Why: Swarm clustering will ot work on Illumina names and
# requires abundance values.
# author: Peter Thorpe and Leighton Pritchard
# September 2016. The James Hutton Insitute, Dundee, UK.

# os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

################################################################
# functions


def dereplicate(seq_dict, abundance_dict, seq_record):
    """function to dereplicate the seq.
    Passed two dictionaries. both dict keys are the coded
    seq_ID e.g. "seqID%d" % (count). seq_dict val is
    the seq_record.seq. abundance_dict val is the abundance"""
    for seq_ID, sequence in seq_dict.items():
        if str(seq_record.seq) == sequence:
            # not a new seq, add one to its abundance
            # of the current key!!
            abundance_dict[seq_ID] += 1
            # it has matched. Break the loop
            return seq_dict, abundance_dict
    # because the loop hasnt broken. This is new:
    # new, create a new dict entry
    seq_dict[seq_record.id] = seq_record.seq
    # set the abundance to 1. First observation
    abundance_dict[seq_record.id] = 1
    return seq_dict, abundance_dict

def return_fasta(seq_dict, abundance_dict, fasta_out):
    """function to convert the two dictionaries into
    a fasta file. seq_dict contains a coded ID e.g.
    "seqID%d" % (count) and the sequence as the val.

    abundance_dict contains the coded ID, as the example
    above "seqID%d" % (count), and the abundance of that
    exact seq
    """
    for seq_ID, sequence in seq_dict.items():
        
        abundance = abundance_dict[seq_ID]
        title = "%s_%d" %(seq_ID, abundance)
        seq_record = SeqRecord(Seq(str(sequence),
                   generic_dna),
                   id=title, name="",
                   description="")
        SeqIO.write(seq_record, fasta_out, "fasta")

def write_out_database(name_out, name):
    """function to write a lists (name) to a file
    name_out (which has already been opended for wrting"""
    #file to keep track of the original names if we need them
    for i in name:
       name_out.write(i)
    name_out.close()
    

def parse_fasta_file(filename, database_out, out, barcode):
    """this function uses biopython to parse the fasta file
    and call other function in order to dereplicate the seq,
    and therefore count the abundance of that seq."""
    fasta_out = open(out, 'w')
    f_in = open(fasta, "r")
    #file to keep track of the original names if we need them
    name_out = open(database_out, "w")

    name = []
    count = 0
    
    seq_dict = dict()
    abundance_dict = dict()
    for seq_record in SeqIO.parse(filename, "fasta"):            
        count = count+1
        if barcode:
            old_to_new_name = "%s\t%sseqID%d\n" % (seq_record.id,
                                                   barcode,
                                                   count)
        else:
            old_to_new_name = "%s\tseqID%d\n" % (seq_record.id,
                                                 count)
        name.append(old_to_new_name)
        #remove the read prefix to get to uniq names
        # underscroe _1 implies an ubundance of 1 for swarm
        if barcode:
            seq_record.id = "%sseqID%d" % (barcode, count)
        else:
            seq_record.id = "seqID%d" % (count)
        seq_record.description = ""
        if count == 1:
            # new, create a new dict entry
            seq_dict[seq_record.id] = seq_record.seq
            # set the abundance to 1. First observation
            abundance_dict[seq_record.id] = 1
            continue
        seq_dict, abundance_dict = dereplicate(seq_dict,
                                               abundance_dict,
                                               seq_record)
    # write out the fasta results
    return_fasta(seq_dict, abundance_dict, fasta_out)
    # write out the coded database 
    write_out_database(name_out, name)
    fasta_out.close()
    f_in.close()
    

#######################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python deduplicate_rename.py -f seq.fasta -d database.out -o out.fasta

script to reformt fasta names.

Name will be reformatted based on SeqID1 - the number is purely
the count of which it is encountened in the file. The abundance is added
as _3, for example.

Barcode seq, if given at the command line will be added to the beginning
of the seq. e.g. ATGTASeqID1_3

and append the abundance of that exact sequence at the end of the
name. For Swarm clustering. 

requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered")

parser.add_option("-d", "--databse", dest="databse", default=None,
                  help="outfile to keep track of old and new ID names",
                  metavar="FILE")
parser.add_option("-b", "--barcode", dest="barcode", default=False,
                  help="barcode seq to add to names on the output")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filename")


(options, args) = parser.parse_args()

fasta = options.fasta
databse = options.databse
barcode = options.barcode
out = options.out

# run the program
parse_fasta_file(fasta, databse, out, barcode)

