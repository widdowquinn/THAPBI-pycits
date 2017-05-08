#!/usr/bin/env python3
#
# santi_pipeline.py
#
# Code for script to identify OTUs from metabarcoding reads.
#
# This is an almost direct translation of a pipeline written by Santiago
# Garcia (Forest Research), to generate OTU clusters from metabarcoding ITS
# reads in Phytophthora.
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe

import logging
import logging.handlers
import multiprocessing
import os
import sys
import time
import traceback

from argparse import ArgumentParser
from biom import load_table

from pycits import tools, fastqc, seq_crumbs, blast, qiime
from pycits import muscle, __version__


# Process command-line arguments
def parse_cmdline():
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="pycits.py")
    parser.add_argument("-p", "--prefix", dest="prefix",
                        action="store", default=None,
                        help="Paired readfiles prefix")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Path to directory containing input reads")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default="script_output",
                        help="Path to directory to write output")
    parser.add_argument("-r", "--reference", dest="reference_fasta",
                        action="store", default=None,
                        help="Path to reference sequence FASTA file")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Report verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-t", "--threads", dest="threads",
                        action="store", type=int,
                        default=multiprocessing.cpu_count(),
                        help="Number of threads to use (default: all)")
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


# Make logger
def construct_logger(args):
    """Returns a logger, configured as required from cmd-line arguments."""
    # Instantiate logger
    logger = logging.getLogger('pycits.py %s: %s' % (time.asctime(),
                                                     __version__))
    logger.setLevel(logging.DEBUG)  # default logger level

    # Set default stream to STDERR
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    # Do we need verbosity in STDERR?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)

    # Was a logfile specified? If so, add it as a stream
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
        except FileNotFoundError:
            logger.error("Could not open %s for logging (exiting)",
                         args.logfile)
            return 1
        else:
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            # logfile is always verbose
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)

    return logger


def report_arguments(logger, cmdline, args):
    """Report usage information to the logger."""
    logger.info("Command-line: %s", cmdline)
    logger.info(args)
    logger.info("Python version: %s", sys.version)
    logger.info("Python executable: %s", sys.executable)
    logger.info("Starting pipeline: %s", time.asctime())


def create_output_directory(outdirname, logger):
    """Create the output directory if it doesn't exist."""
    try:
        os.makedirs(outdirname)
    except OSError:
        logger.error("Error creating directory %s (exiting)",
                     outdirname)
        logger.info(last_exception())
        return 1
    else:
        if not os.path.isdir(outdirname):
            logger.error("%s is not a directory (exiting)", outdirname)
            return 1
    return 0


def get_filenames(args, logger):
    """Return input filenames, or 1 if any input files contain whitespace.

    Whitespace causes problems with some of the third-party applications, so
    we want to avoid it.
    """
    infilenames = sorted([os.path.join(args.indirname, fname) for
                          fname in os.listdir(args.indirname) if
                          fname.startswith(args.prefix)])
    logger.info("Input files: %s", infilenames)
    for fname in infilenames:
        if ' ' in os.path.relpath(fname):
            logger.error("Relative path to file or directory " +
                         "'%s' contains whitespace (exiting)", fname)
            return 1
    return infilenames


def create_interfaces(args, logger):
    """Check third-party executable dependencies exist and instantiate.

    The instantiated interfaces are returned.
    """
    # Check for presence of third-party tools, by instantiating interfaces
    logger.info("Checking third-party packages:")
    logger.info("\tFastQC... (%s)", args.fastqc)
    fastqc_qc = fastqc.FastQC(args.fastqc)
    logger.info("\ttrim_quality... (%s)", args.trim_quality)
    trim_quality = seq_crumbs.Trim_Quality(args.trim_quality, logger)
    logger.info("\tjoin_paired_ends.py... (%s)", args.join_paired_ends)
    jpe = qiime.Join_Paired_Ends(args.join_paired_ends, logger)
    logger.info("\tconvert_format... (%s)", args.convert_format)
    convert_format = seq_crumbs.Convert_Format(args.convert_format, logger)
    logger.info("\tblastclust... (%s)", args.blastclust)
    blastclust = blast.Blastclust(args.blastclust)
    logger.info("\tmuscle... (%s)", args.muscle)
    muscle_aln = muscle.Muscle(args.muscle)
    logger.info("\tpick_otus.py... (%s)", args.pick_otus)
    pick_otus = qiime.Pick_Otus(args.pick_otus, logger)
    logger.info("\tpick_closed_reference_otus.py... (%s)",
                args.pick_closed_reference_otus)
    pcro = qiime.Pick_Closed_Ref_Otus(args.pick_closed_reference_otus, logger)
    return (fastqc_qc, trim_quality, jpe, convert_format, blastclust,
            muscle_aln, pick_otus, pcro)


def trim_reads(trim_quality, infilenames, args, logger):
    """Trim input reads for quality and return."""
    logger.info("Trim reads by quality")
    trimmed_fnames = [trim_quality.run(fname, args.outdirname) for
                      fname in infilenames]
    logger.info("Trimmed FASTQ files:")
    return trimmed_fnames


def join_reads(jpe, trimmed_fnames, args, logger):
    logger.info("Join trimmed, paired-end reads")
    joined_reads = jpe.run(trimmed_fnames, args.outdirname)
    logger.info("Joined reads:\t%s", joined_reads)
    return joined_reads


def fastq_to_fasta(convert_format, infname, outdirname, logger):
    """Convert FASTQ sequence to FASTA."""
    logger.info("Creating FASTA file from %s", infname)
    fasta = convert_format.run(infname, outdirname)
    logger.info("Converted to FASTA:\t%s", fasta)
    return fasta


def trim_sequences(infname, outfname, logger):
    """Trim sequences by 20nt on left and right."""
    logger.info("Trimming sequences")
    write_count = tools.trim_seq(infname, outfname)
    logger.info("Trimmed, joined FASTA (%d sequences):\t%s",
                write_count, outfname)


def blastclust_cluster(blastclust, infname, outdirname, threads, logger):
    """Cluster sequences with BLASTCLUST."""
    logger.info("Clustering OTUs with BLASTCLUST")
    bc_result = blastclust.run(infname, outdirname, threads)
    logger.info("Clustering joined, trimmed sequences with " +
                "BLASTCLUST:\n\t%s", bc_result.outfilename)
    return bc_result


def muscle_align(muscle_aln, infiles, outdir, logger):
    """Align files with MUSCLE and place alignments in outdir."""
    logger.info("Aligning BLASTCLUST OTU sequences with MUSCLE")
    for fname in infiles:
        muscle_result = muscle_aln.run(os.path.join(fname))
    logger.info("Aligned BLASTCLUST OTU sequences written to:\n\t%s",
                outdir)
    return muscle_result


def pick_denovo_otus(pick_otus, infname, reference, outdirname, logger):
    """Pick denovo OTUs from input sequences with QIIME."""
    logger.info("Picking UCLUST OTUs with QIIME")
    qiime_uclustdir = pick_otus.run(infname, reference, outdirname)
    logger.info("OTUs picked by QIIME with UCLUST written to:\n\t%s",
                qiime_uclustdir)
    return qiime_uclustdir


def pick_closedref_otus(pcro, infname, reference, outdirname, logger):
    """Pick closed-reference OTUs with QIIME."""
    logger.info("Picking closed-reference OTUs with QIIME")
    qiime_pcrodir = pcro.run(infname, reference, outdirname)
    logger.info("OTUs picked by QIIME (closed-reference) " +
                "written to:\n\t%s", qiime_pcrodir)
    return qiime_pcrodir


def biom_to_tsv(biomfname, logger):
    """Converts BIOM file to TSV and places output in same directory."""
    logger.info("Converting BIOM output to tabular (TSV) format")
    biom_table = load_table(biomfname)
    tsvfname = os.path.splitext(biomfname)[0] + ".tsv"
    with open(tsvfname, 'w') as ofh:
        ofh.write(biom_table.to_tsv())
    logger.info("TSV output written to:\n\t%s", tsvfname)
    return tsvfname


def run_fastqc(fastqc_qc, infnames, outdirname, logger):
    logger.info("Running FastQC")
    fastqc_outdir = os.path.join(outdirname, "FastQC")
    results = []
    for infname in infnames:
        logger.info(infname)
        qc_result = fastqc_qc.run(infname, fastqc_outdir)
        logger.info("Writing to:\n\t%s\n\t%s",
                    qc_result.htmlfile, qc_result.zipfile)
        results.append(qc_result)
    return results


# Main process for script
def run_pycits_main(namespace=None):
    """Main process for run-pycits.py script"""
    # Parse command-line if needed
    if namespace is None:
        args = parse_cmdline()
    else:
        args = namespace

    # Set up logging
    logger = construct_logger(args)

    # Report input arguments
    report_arguments(logger, ' '.join(sys.argv), args)

    # Have we got an input directory, reference set and prefix? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        return 1
    logger.info("Input directory: %s", args.indirname)

    if args.prefix is None:
        logger.error("No read file prefix given (exiting)")
        return 1
    logger.info("Read file prefix: %s", args.prefix)

    if args.reference_fasta is None:
        logger.error("No reference FASTA file given (exiting)")
        return 1
    logger.info("Reference FASTA file: %s", args.reference_fasta)

    # Have we got an output directory and prefix? If not, create it.
    # create_output_directory() returns 1 if creation fails.
    logger.info("Creating directory %s", args.outdirname)
    if create_output_directory(args.outdirname, logger):
        return 1

    # Check for the presence of space characters in any of the input filenames
    # If we have any, abort here and now.
    infilenames = get_filenames(args, logger)
    if infilenames is None:
        return 1

    # Check for presence of third-party tools, by instantiating interfaces
    (fastqc_qc, trim_quality, jpe, convert_format, blastclust,
     muscle_aln, pick_otus, pcro) = create_interfaces(args, logger)

    # How many threads are we using?
    args.threads = min(args.threads, multiprocessing.cpu_count())
    logger.info("Using %d threads/CPUs where available", args.threads)

    # Trim reads on quality - forward and reverse reads
    trimmed_fnames = trim_reads(trim_quality, infilenames, args, logger)

    # Join the trimmed, paired-end reads together
    joined_reads = join_reads(jpe, trimmed_fnames, args, logger)

    # Create a FASTA file equivalent to the joined FASTQ reads
    joined_fasta = fastq_to_fasta(convert_format, joined_reads,
                                  args.outdirname, logger)

    # Trim sequences by 20nt on left and right
    trimmed_joined_fasta = os.path.splitext(joined_fasta)[0] + '_trimmed.fasta'
    write_count = trim_sequences(joined_fasta, trimmed_joined_fasta, logger)

    # Cluster OTUs with BLASTCLUST
    bc_result = blastclust_cluster(blastclust, trimmed_joined_fasta,
                                   args.outdirname, args.threads, logger)

    # Convert BLASTCLUST output to FASTA sequence files
    logger.info("Generating FASTA from BLASTCLUST output")
    bc_outdir, bc_files = tools.blastclust_to_fasta(bc_result.outfilename,
                                                    trimmed_joined_fasta,
                                                    args.outdirname)
    logger.info("FASTA sequences for BLASTCLUST OTUs written to:\n\t%s",
                bc_outdir)

    # Align the BLASTCLUST OTUs with MUSCLE
    muscle_result = muscle_align(muscle_aln, bc_files, bc_outdir, logger)

    # Pick de novo OTUs with QIIME
    qiime_uclustdir = pick_denovo_otus(pick_otus, trimmed_joined_fasta,
                                       args.reference_fasta, args.outdirname,
                                       logger)

    # Pick closed-reference OTUs with QIIME
    qiime_pcrodir = pick_closedref_otus(pcro, trimmed_joined_fasta,
                                        args.reference_fasta, args.outdirname,
                                        logger)

    # Convert BIOM output to TSV
    biomfname = os.path.join(qiime_pcrodir, "otu_table.biom")
    tsvfname = biom_to_tsv(biomfname, logger)

    # Run FastQC on the read files
    qc_results = run_fastqc(fastqc_qc,
                            infilenames + trimmed_fnames + [joined_reads],
                            args.outdirname, logger)

    # Announce end of pipeline
    logger.info("Pipeline complete: %s", time.asctime())
    return 0
