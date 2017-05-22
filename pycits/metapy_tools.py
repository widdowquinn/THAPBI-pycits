#!/usr/bin/env python3
#
# metapy_tools.py
#
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import gzip
from .tools import convert_fq_to_fa
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
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
    return os.path.splitext(infile)[0]


def test_reads_exist_and_suffix(reads):
    """function to test if the reads exist and if there
    are compressed or not.
    Take in a read file.
    This is required, if a run is cancelled then re-run
    with the same command and the reads have been decompresssed,
    it will fail.
    return the proper name the reads will have"""
    if not os.path.isfile(reads):
        # there is something wrong ...
        # this file doesnt exist. Maybe already decompressed
        # but not what the user is asking for
        if os.path.isfile(os.path.splitext(reads)[0]):
            # the reads have already been decompressed
            return os.path.splitext(reads)[0]
    if os.path.isfile(reads):  # file is real
        if reads.endswith(".gz"):  # decompress them
            new_name = decompress(reads)
            return new_name
        else:
            return reads
    else:
        error = "\nERROR: %s   FILE DOES NOT EXIT" % reads
        sys.exit(error + ". Check your input file path and name\n")


def make_folder(folder, WORKING_DIR, exist_ok=True):
    """function to make a folder with desired name"""
    dest_dir = os.path.join(WORKING_DIR, folder)
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists " +
               "I will write over what is in there!!")
    return dest_dir


def compress(infile):
    """function to compress reads, make them .gz"""
    cmd = ' '.join(["gzip", "-f", infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)


def covert_chop_read(infile, LEFT_TRIM, RIGHT_TRIM):
    """function to reduce repetive code:
    Take in an assembled fq file, either PEAR or FLASH
    outfile. Converts this to Fasta, then chops the seq
    at LEFT and RIGHT
    write out: infile + '.bio.fasta'
    infile + '.bio.chopped.fasta """
    convert_fq_to_fa(infile,
                     infile + ".bio.fasta")
    # need to trim the left and right assembled seq so they
    # cluster with the database.
    # use: trim_seq() from tools.
    # trim_seq(infname, outfname, lclip=53, rclip=0, minlen=100)
    # this function is now added to this script
    metapy_trim_seq(infile + ".bio.fasta",
                    infile + ".bio.chopped.fasta",
                    LEFT_TRIM, RIGHT_TRIM)

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


def database_checker(filename, outfile="temp.fasta"):
    """this function re-write a file as a fasta file.
    I am using this to test if the DB is ill formatted.
    Biopython should detect this if it is badly formatted.
    No duplicate name or sequences are allowed.
    These are detected by sets.
    Now checks for illegal characters. Only allowed:
    ATCG

    Returns: error, problem seqrecord

    or

    return "ok", "ok" if passed checks"""
    f = open(outfile, 'w')
    name_set = set([])
    seq_set = set([])
    nt_set = set("ATCG")
    if not os.path.isfile(filename):
        error = "\nDATABASE: %s does NOT EXIST. Check " % filename
        sys.exit(error + " the correct name  and PATH\n")
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = str(seq_record.seq)
        # check for illegal characters
        if set(seq.upper()).difference(nt_set):
            illegal = set(str(seq_record.seq)).difference(nt_set)
            error_str = "Illegal nucleotide found %s. " % illegal
            err = error_str + "Only %s are allowed. " % (str(nt_set))
            return err + " Problem in: ", seq_record
        # check for duplicate names. Should be unique
        if seq_record.id not in name_set:
            name_set.add(seq_record.id)
        else:
            return "Duplicate names found", seq_record
        # check duplicate sequences. should be unique
        if str(seq_record.seq) not in seq_set:
            seq_set.add(str(seq_record.seq))
        else:
            return "Duplicate sequence found", seq_record
        try:
            SeqIO.write(seq_record, f, "fasta")
        except ValueError:
            return "Ill formatted fasta file", seq_record
    f.close()
    return "ok", "ok"
