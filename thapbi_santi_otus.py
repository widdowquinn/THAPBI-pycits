#!/usr/bin/env python
#
# thapbi_santi_otus.py
#
# Script to identify OTUs from metabarcoding reads.
#
# This is an almost direct translation of a pipeline written by Santiago
# Garcia, to generate OTU clusters from metabarcoding ITS reads in
# Phytophthora.
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import logging
import logging.handlers
import os
import sys
import traceback

from argparse import ArgumentParser

from thapbi_santi import fastqc, seq_crumbs, ea_utils, blast, qiime

# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="thapbi_santi_otus.py")
    parser.add_argument("-p", "--prefix", dest="prefix",
                        action="store", default=None,
                        help="Paired readfiles prefix")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Path to directory containing input reads")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default="script_output",
                        help="Path to directory to write output")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Report verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--fastqc", dest="fastqc",
                        action="store", default="fastqc",
                        help="Path to FastQC executable")
    parser.add_argument("--trim_quality", dest="trim_quality",
                        action="store", default="trim_quality",
                        help="Path to seq_crumbs trim_quality script")
    parser.add_argument("--join_paired_ends", dest="join_paired_ends",
                        action="store", default="join_paired_ends.py",
                        help="Path to ea-utils join_paired_ends.py script")
    parser.add_argument("--convert_format", dest="convert_format",
                        action="store", default="convert_format",
                        help="Path to seq_crumbs convert_format script")
    parser.add_argument("--blastclust", dest="blastclust",
                        action="store", default="blastclust",
                        help="Path to blastclust")
    parser.add_argument("--muscle", dest="muscle",
                        action="store", default="muscle",
                        help="Path to MUSCLE")
    parser.add_argument("--pick_otus", dest="pick_otus",
                        action="store", default="pick_otus.py",
                        help="Path to QIIME pick_otus.py script")
    parser.add_argument("--pick_closed_reference_otus",
                        dest="pick_closed_reference_otus",
                        action="store",
                        default="pick_closed_reference_otus.py",
                        help="Path to QIIME pick_closed_reference_otus.py " +
                        "script")
    return parser.parse_args()

# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('thapbi_santi_otus.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)

    # Have we got an input directory and prefix? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.indirname)
    if args.prefix is None:
        logger.error("No read file prefix given (exiting)")
        sys.exit(1)

    # Have we got an output directory and prefix? If not, create it.
    if not os.path.exists(args.outdirname):
        logger.warning("Output directory %s does not exist - creating it" %
                       args.outdirname)
        os.makedirs(args.outdirname)

    # Check for the presence of space characters in any of the input filenames
    # If we have any, abort here and now.
    infilenames = sorted([os.path.join(args.indirname, fname) for 
                          fname in os.listdir(args.indirname) if
                          fname.startswith(args.prefix)])
    logger.info("Input files: %s" % infilenames)
    for fname in infilenames:
        if ' ' in  os.path.abspath(fname):
            logger.error("File or directory '%s' contains whitespace" % fname)
            logger.error("(exiting)")
            sys.exit(1)

    # Check for presence of third-party tools, by instantiating interfaces
    logger.info("Checking third-party packages:")
    logger.info("\tFastQC... (%s)" % args.fastqc)
    fastQC = fastqc.FastQC(args.fastqc, logger)
    logger.info("\ttrim_quality... (%s)" % args.trim_quality)
    trim_quality = seq_crumbs.Trim_Quality(args.trim_quality, logger)
    logger.info("\tjoin_paired_ends.py... (%s)" % args.join_paired_ends)
    jpe = qiime.Join_Paired_Ends(args.join_paired_ends, logger)
    logger.info("\tconvert_format... (%s)" % args.convert_format)
    convert_format = seq_crumbs.Convert_Format(args.convert_format, logger)
    logger.info("\tblastclust... (%s)" % args.blastclust)
    blastclust = blast.Blastclust(args.blastclust, logger)
    logger.info("\tpick_otus.py... (%s)" % args.pick_otus)
    pick_otus = qiime.Pick_Otus(args.pick_otus, logger)
    logger.info("\tpick_closed_reference_otus.py... (%s)" %
                args.pick_closed_reference_otus)
    pcro = qiime.Pick_Closed_Ref_Otus(args.pick_closed_reference_otus, logger)
    
    # Trim reads on quality - forward and reverse reads
    logger.info("Trim reads by quality")
    try:
        trimmed_fnames = [trim_quality.run(fname, args.outdirname) for
                          fname in infilenames]
        logger.info("Trimmed FASTQ files:")
        logger.info("\t%s" % trimmed_fnames)
    except:
        logger.error("Error running trim_quality (exiting)")
        logger.error(last_exception())
        sys.exit(1)
    
    # Join the trimmed, paired-end reads together
    logger.info("Join trimmed, paired-end reads")
    try:
        joined_reads = jpe.run(trimmed_fnames, args.outdirname)
        logger.info("Joined reads:")
        logger.info("\t%s" % joined_reads)
    except:
        logger.error("Error joining reads (exiting)")
        logger.error(last_exception())
        sys.exit(1)

    # Create a FASTA file equivalent to the joined FASTQ reads
    logger.info("Creating FASTA file")
    try:
        joined_fasta = convert_format.run(joined_reads, args.outdirname)
        logger.info("Converted to FASTA:")
        logger.info("\t%s" % joined_fasta)
    except:
        logger.error("Error converting to FASTA (exiting)")
        logger.error(last_exception())
        sys.exit(1)
