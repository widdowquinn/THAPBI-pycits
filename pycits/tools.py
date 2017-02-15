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
import pysam

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
    abundance = defaultdict(int)
    hash_to_seq = defaultdict(str)
    hash_to_name = defaultdict(str)
    new_to_old = defaultdict(list)
    # compile hashed sequence data - the loop is necessary because we need
    # to run through completely once to obtain abundance info. If memory
    # to store this becomes an issue, we might want to try an alternative
    # approach
    for seq in seqlist:
        seqhash = hashlib.md5(str(seq.seq).encode()).hexdigest()
        abundance[seqhash] += 1
        hash_to_seq[seqhash] = seq.seq
        hash_to_name[seqhash] = seq.id
        new_to_old[seqhash].append(seq.id)
    # generate list of sequences, named for the hash, including abundance data,
    # and return with ID:hash lookup
    seqlist = list()
    for name, abundance_val in abundance.items():
        seqname = "{0}_{1}".format(name, abundance_val)
        seqlist.append(SeqRecord(id=seqname, description="",
                                 seq=hash_to_seq[name]))
    return(seqlist, new_to_old, hash_to_seq)


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
    seqlist, new_to_old, hash_to_seq = deduplicate_and_rename(seqlist)
    for key, vals in new_to_old.items():
        out_data = "%s\t%s\n" % (key, "\t".join(vals))
        name_out.write(out_data)
    for seq_record in (seqlist):
        SeqIO.write(seq_record, fasta_out, "fasta")
    # close the open files
    fasta_out.close()
    name_out.close()


def check_OTU_db_abundance_val(OTUfasta):
    """this is a function to check that the OTU database has _1
    at the end of the title, this is required for swamr clustering
    takes in a fasta and if required, writes out a new one"""
    # we will iterate through first to see if we need to rewrite it.
    bad_count = 0
    for seq_record in SeqIO.parse(OTUfasta, "fasta"):
        if not seq_record.id.endswith("_1"):
            # print ("OTU db %s missing abundance value" % seq_record.id)
            bad_count = bad_count + 1
    if bad_count > 0:
        # only if we need to write to a new file, we will.
        outfile = OTUfasta.split(".fa")[0] + "abundance.fasta"
        f_out = open(outfile, 'w')
        for seq_record in SeqIO.parse(OTUfasta, "fasta"):
            seq_record.id = seq_record.id + "_1"
            seq_record.description = ""
            SeqIO.write(seq_record, f_out, "fasta")
        f_out.close()
        return outfile
    else:
        return OTUfasta


def filter_sam_file(in_sam, outfile):
    """function for pysam to filter the sam file to obtain matches
    that only match for its length
    # AS:i:0 are perfect matches
    # samfile format, good apge:
    # http://www.metagenomics.wiki/tools/samtools/bam-sam-file-format
    """
    # counter to track number of seq of interest
    no_missmtch = open(outfile, "w")
    samfile = pysam.AlignmentFile(in_sam)
    cig_list = []
    matches = []
    for read in samfile:
        # this apparently are pure matches. But not the same as AS:i:0(?)
        if (len(read.cigar)) == 1:
            out_info = "%s\t%s\t%s\n" % (read.reference_name,
                                         read.query_name,
                                         read.cigar)
            no_missmtch.write(out_info)
            matches.append(out_info)
        cig_list.append(read.cigar)
    no_missmtch.close()
    return cig_list, matches


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


# the following four function are to rename the clusters back to their
# original names
def get_names_from_Seq_db(seq_db):
    """function to get a list of name in the seq db"""
    names = []
    names_abudance_removed = []
    for seq_record in SeqIO.parse(seq_db, "fasta"):
        if seq_record.id.endswith("_1"):
            names.append(seq_record.id)
        else:
            names_abudance_removed.append(seq_record.id)
            names.append(seq_record.id + "_1")
    return names, names_abudance_removed


def coded_name_to_species(database_file):
    """functiong takes the already generated tab separated
    database of coded name to species file. Returns a dic
    of coded_name to species"""
    with open(database_file) as file:
        data = file.read().split("\n")
    coded_name_to_species_dict = dict()
    for line in data:
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith("#"):
            continue
        data = line.split("\t")
        coded_name = data[0]
        species = data[1:]
        coded_name_to_species_dict[coded_name.rstrip()] = species
    return coded_name_to_species_dict


def return_real_line(line):
    """function to return only true lines,
    no comments, or blank lines."""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):  # dont want comment lines
        return False
    if "\t" in line:
        cluster_line = line.rstrip("\n").split("\t")
    else:
        # different clustering program?
        cluster_line = line.rstrip("\n").split()
    return cluster_line


def parse_tab_file_get_clusters(in_file, seq_db, database, out_file):
    """script to open up a tab or space separeted clustering
    output and rename according to the name in the database file.
    Abundance is also appended to the name"""
    # call the function to get the dictionary
    # populated with the database
    names, names_abudance_removed = get_names_from_Seq_db(seq_db)
    coded_name_to_species_dict = coded_name_to_species(database)
    cluster_file = open(in_file, "r")
    summary_out_file = open(out_file, "w")
    count = int(0)
    for line in cluster_file:
        output_str = ""
        count += 1
        # get func to return real line only
        cluster_line = return_real_line(line)
        if not cluster_line:
            continue
        for member in cluster_line:
            if member in names or member in names_abudance_removed:
                if member.endswith("_1"):
                    member = ("_").join(member.split("_")[:-1])
                cluster_summary = "%s_abundance=1\t" % (member)
                output_str = output_str + cluster_summary
                continue
            try:
                # data would be: 2f1454d16278fda2d44f26ebf7a0ed05_310
                # _abundance value
                split_name = member.split("_")[:-1]
                coded_name = ("_").join(split_name)
                if coded_name in names:
                    species = coded_name
                    abundance = "1"
                else:
                    species = coded_name_to_species_dict[coded_name]
                    abundance = member.split("_")[-1]
                    species = "\t".join("%s_abundance=%s"
                                        % (i, abundance) for i in species)
                # print (abundance)
            except:
                KeyError
                print ("we have an error at %s " % member)
                print ("""something went wrong with decoding the names maybe
                your names were not separated in your database? If so, the
                block of code above this staement need adjusting accordingly.
                """)
                sys.exit()
            # add the info to a str, we will write at the
            # end of the cluster line
            cluster_summary = "%s\t" % (species)
            output_str = output_str + cluster_summary
        summary_out_file.write(output_str+"\n")
    # close the files
    cluster_file.close()
    summary_out_file.close()
