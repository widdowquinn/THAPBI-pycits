#!/usr/bin/env python

import os

from Bio import SeqIO
from subprocess import check_output

def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        exefile = check_output(["which", filename]).strip()
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)

# Function replacing Santi's trim_longitudes.py script
# Defaults are equivalent to the default settings in Santi's script, excepting
# the filenames, which we now provide directly.
def trim_seq(infname, outfname, lclip=21, rclip=20, minlen=100):
    """Trims all FASTA sequences in infname by lclip and rclip on the
    'left'- and 'right'-hand ends, respectively and writes the results to
    outfname. Any clipped sequence with length < minlen is not written.
    """
    with open(infname, 'r') as fh:
        # Use generators to save memory
        s_trimmed = (s[lclip:-rclip] for s in SeqIO.parse(fh, 'fasta'))
        return SeqIO.write((s for s in s_trimmed if len(s) >= minlen),
                           outfname, 'fasta')

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
            
