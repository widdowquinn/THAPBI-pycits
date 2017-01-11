#!/usr/bin/env python
#
# tools.py
#
# Reimplementations of Santi's trim_longitudes.py and blastclust_lst2fasta.py
# scripts/functions, and other miscellaneous functions
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os

import gzip
from subprocess import check_output, CalledProcessError
import hashlib
import sys

from collections import defaultdict
from optparse import OptionParser

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class NotExecutableError(Exception):
    """Exception raised when expected executable is not executable"""
    def __init__(self, message):
        self.message = message


def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        try:
            exefile = check_output(["which", filename]).strip()
        except CalledProcessError:
            raise NotExecutableError("{0} does not exist".format(filename))
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)


# Defaults  settings are for the current ITS1 phy
# primers and the development ITS1 Phy database. DON'T BLINDLY USE THESE
# The filenames, we now provide directly.
def trim_seq(infname, outfname, lclip=53, rclip=0, minlen=100):
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


def deduplicate_and_rename(seqlist):
    """Removes duplicates from the passed SeqRecords, and replaces sequence
    IDs with the md5 hash of the sequence, suffixed with the number of
    sequences that were replaced. Returns a tuple (sequences, names_old_to_new)
    where names_old_to_new is a list of tuples of (old_name, hash) pairs.

    - sequences   Iterable of Bio.SeqRecord objects
    """
    names_old_to_new = list()
    abundance = defaultdict(int)
    hash_to_seq = defaultdict(str)
    hash_to_name = defaultdict(str)
    # compile hashed sequence data - the loop is necessary because we need
    # to run through completely once to obtain abundance info. If memory
    # to store this becomes an issue, we might want to try an alternative
    # approach
    for seq in seqlist:
        seqhash = hashlib.md5(str(seq.seq).encode()).hexdigest()
        abundance[seqhash] += 1
        hash_to_seq[seqhash] = seq.seq
        hash_to_name[seqhash] = seq.id
        database_output = ("{0}\t{1}\n".format(seq.id, seqhash))
        names_old_to_new.append(database_output)
    # generate list of sequences, named for the hash, including abundance data,
    # and return with ID:hash lookup
    seqlist = list()
    for name, abundance_val in abundance.items():
        seqname = "{0}_{1}".format(name, abundance_val)
        seqlist.append(SeqRecord(id=seqname, description="",
                                 seq=hash_to_seq[name]))
    return(seqlist, names_old_to_new, hash_to_seq)


def dereplicate_name(fasta, database_out, out):
    """function to dereplicate the seq. Thus generating an
    abundance of that seq. Rename the seq to something
    swarm will work with.
    - fasta     - fasta file of assembled seq
    - database_out   - ouput old name to new file.
    -out        - outfile
    """
    # open files to write to
    fasta_out = open(out, 'w')
    name_out = open(database_out, "w")
    # convert the file to a list of Seqrecord objects
    seqlist = list(SeqIO.parse(fasta, 'fasta'))
    seqlist, names_old_to_new, hash_to_seq = deduplicate_and_rename(seqlist)
    for i in (names_old_to_new):
        name_out.write(i)
    for seq_record in (seqlist):
        SeqIO.write(seq_record, fasta_out, "fasta")
    # close the open files
    fasta_out.close()
    name_out.close()


# Function replacing Santi's blastclust_lst2fasta.py script
def blastclust_to_fasta(infname, seqfname, outdir):
    """Converts input BLASTCLUST output list to a subdirectory of FASTA files,
    each of which contains all sequences from a single cluster. The sequences
    matching the IDs listed in BLASTCLUST output should all be in the file
    seqfname.
    """
    outdirname = os.path.join(outdir, "blastclust_OTUs")
    if not os.path.exists(outdirname):
        os.makedirs(outdirname)
    seqdict = SeqIO.index(seqfname, 'fasta')
    with open(infname, 'r') as fh:
        otu_id = 0
        for line in fh:
            otu_id += 1
            outfname = os.path.join(outdirname,
                                    "blastclust_OTU_%06d.fasta" % otu_id)
            SeqIO.write((seqdict[key] for key in line.split()),
                        outfname, 'fasta')
    return outdirname
