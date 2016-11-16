#!/usr/bin/env python
# THIS IS a tester for module dev only.
#
# TODO - put this functionality under a proper test framework (e.g. nosetests
#        or equivalent)


import errno
import logging
import logging.handlers
import multiprocessing
import os
import sys
import time
import traceback

from argparse import ArgumentParser

from pycits import tools, trimmomatic, pear, error_correction,\
    clean_up, swarm, seq_crumbs, deduplicate, \
    rename_clusters_new_to_old, chop_seq

# setting up some test variables
# LP NOTE: globals should be capitalised
THREADS = "4"  # threads to use

# test read set
LEFT_READS = "./data/reads/DNAMIX_S95_L001_R1_001.fastq"
RIGHT_READS = "./data/reads/DNAMIX_S95_L001_R2_001.fastq"
READ_PREFIX = os.path.split(LEFT_READS)[-1].split("_R")[0]
assert(READ_PREFIX == os.path.split(RIGHT_READS)[-1].split("_R")[0])

# left_trim 53 - this should be default??
LEFT_TRIM = 53
RIGHT_TRIM = 0

# Third-party applications
SPADES_EXE = "/home/pt40963/scratch/Downloads/SPAdes-3.9.0-Linux/bin/spades.py"
TRIMMO_JAR = "/home/pt40963/scratch/Downloads/Trimmomatic-0.36/" +\
             "trimmomatic-0.36.jar"
PEAR_EXE = "pear"


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Run as script
if __name__ == '__main__':
    # Set up logging
    logger = logging.getLogger('thapbi_santi_otus.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    try:
        logstream = open("module_testing.txt", 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     args.logfile)
        sys.exit(1)

    # Report input arguments
    logger.info("Command-line: %s" % ' '.join(sys.argv))
    logger.info("Starting testing: %s" % time.asctime())

    # trimmomatic testing.
    logger.info("starting trimmomatic testing")
    trim = trimmomatic.Trimmomatic(TRIMMO_JAR, logger)

    logger.info("Trim reads by quality")
    trim.run(TRIMMO_JAR, LEFT_READS, RIGHT_READS,
             THREADS, "trimmed_reads", 0, logger)

    # error correction testing
    error_correct = error_correction.Error_correction(SPADES_EXE, logger)
    left_trimmed = os.path.join("trimmed_reads",
                                READ_PREFIX +
                                "_paired_R1.fq.gz")
    right_trimmed = os.path.join("trimmed_reads",
                                 READ_PREFIX +
                                 "_paired_R2.fq.gz")
    logger.info("error correction using Bayes hammer")
    error_correct.run(SPADES_EXE, left_trimmed, right_trimmed,
                      "error_correction", THREADS, logger)
    left_corrected = os.path.join("error_correction", "corrected",
                                  READ_PREFIX + "_paired_R1.fq" +
                                  ".00.0_0.cor.fastq.gz")

    right_corrected = os.path.join("error_correction", "corrected",
                                   READ_PREFIX + "_paired_R2.fq" +
                                   ".00.0_0.cor.fastq.gz")

    # PEAR testing - assemble
    assemble = pear.Pear(PEAR_EXE, logger)
    logger.info("assembly using PEAR")
    assemble.run(left_trimmed, right_trimmed,
                 THREADS, "PEAR_assembled", logger)

    # Pear testing with EC reads
    logger.info("assembly using PEAR with error corrected reads")
    assemble.run(left_corrected, right_corrected,
                 THREADS, "PEAR_assembled_EC", logger)

    # cleaning up unwanted files
    logger.info("testing the removal of unwanted file")
    unwanted_list = ["unpaired_R*.fq.gz",
                     "PEAR_assembled/*.discarded*fastq",
                     "PEAR_assembled/*.unassembled*fastq"]

    delete_file = clean_up.Clean_up("rm", logger)
    for i in unwanted_list:
        delete_file.run(i)

    # convert format
    format_change = seq_crumbs.Convert_Format("convert_format",
                                              logger)
    assembled_reads = os.path.join("PEAR_assembled",
                                   READ_PREFIX +
                                   "_paired_PEAR.assembled.fastq")

    format_change.run(assembled_reads, "fasta_converted",
                      logger)

    # deduplicate read:
    #python deduplicate_rename.py
    #-f DNAMIX_S95_L001_paired_PEAR.assembled.fasta
    #-d database.out
    #-o temp.fasta
    logger.info("deduplicating reads")
    dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
    dedup = deduplicate.Deduplicate(dedup_prog, logger)
    fasta_file = os.path.join("fasta_converted",
                              "DNAMIX_S95_L001_paired_PEAR.assembled.fasta")
    database = os.path.join("fasta_converted", "database.out")
    out_dedup = os.path.join("fasta_converted", "temp.fasta")
    dedup.run(dedup_prog, fasta_file, database, out_dedup,
              logger)

    # need to trim the left and right assembled seq so they
    # cluster with the database.
    chop_prog = os.path.join("pycits", "trim_fasta_file.py")

    fasta_file = os.path.join("fasta_converted", "temp.fasta")
    ITS_database = os.path.join("data",
                                "ITS_database_NOT_confirmed_correct_" +
                                "last14bases_removed.fasta")
    fasta_out = os.path.join("fasta_converted",
                             "database_assembled_reads.fasta")

    chopper = chop_seq.Chop_seq(chop_prog, logger)
    #        run(exe_path, fasta, database, left, right,
            #fasta_out, barcode="0", logger=False)
    chopper.run(chop_prog, fasta_file, ITS_database, LEFT_TRIM,
                RIGHT_TRIM, fasta_out, logger)

    # SWARM testing - assemble
    cluster = swarm.Swarm("swarm", logger)
    assembled_fa_reads = os.path.join("fasta_converted",
                                      "temp.fasta")

    logger.info("clustering with Swarm")
    cluster.run(assembled_fa_reads,
                THREADS, 1, "Swarm_cluster_assembled_reads_only", logger)
    cluster.run(fasta_out,
                THREADS, 1, "Swarm_cluster", logger)

    # recode cluster output
    #python pycits/parse_clusters_new_to_old_name.py -i
    # -i swarm_tempd1.out -d database.out -o decoded_clusters.out
    logger.info("renaming swarm output")
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    rename = rename_clusters_new_to_old.Rename(rename_prog, logger)
    infile = os.path.join("Swarm_cluster", "swarm_clustering_d1",
                          "swarm_clustering_d1")
    outfile = os.path.join("Swarm_cluster", "swarm_clustering_d1",
                           "swarm_clustering_d1_renamed.out")
    rename.run(rename_prog, infile, ITS_database,
               database, outfile, logger)

    # need to summarise the clusters now
    # how? scripts already made, are they good enough?
