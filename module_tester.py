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
import subprocess

from argparse import ArgumentParser

from pycits import tools, trimmomatic, pear, error_correction,\
    flash, clean_up, swarm, seq_crumbs#, deduplicate, \
    #rename_clusters_new_to_old, chop_seq

# setting up some test variables
# LP NOTE: globals should be capitalised
THREADS = "4"  # threads to use

# test read set
LEFT_READS = "./data/reads/DNAMIX_S95_L001_R1_001.fastq"
RIGHT_READS = "./data/reads/DNAMIX_S95_L001_R2_001.fastq"
READ_PREFIX = os.path.split(LEFT_READS)[-1].split("_R")[0]
PHREDSCORE = "phred33"
PREFIX = READ_PREFIX
ADAPTERS = os.path.join("adapters", "TruSeq3-PE.fa")
OUTDIR_TRIM = "trimmed_reads"
OUTFILES = [os.path.join(OUTDIR_TRIM, PREFIX + suffix) for suffix in
           ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
            "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]


assert(READ_PREFIX == os.path.split(RIGHT_READS)[-1].split("_R")[0])

#####################################################################
# make a folder_for_all_the_out_files
folder_list = ["trimmed_reads", "PEAR_assembled", "PEAR_assembled_EC",
   	       "Swarm_cluster_assembled_reads_only", "Swarm_cluster",
               "flash_assembled", "flash_assembled_EC"]
file_name = 'test.txt'
working_dir = os.getcwd()
for i in folder_list:
	dest_dir = os.path.join(working_dir, i)
	try:
    		os.makedirs(dest_dir)
	except OSError:
		print ("""folder already exists, 
                      I will write over what is in there!!""")

########################################################################
# left_trim 53 - this should be default??
LEFT_TRIM = 53
RIGHT_TRIM = 0


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))

########################################################################
# Run as script
if __name__ == '__main__':
    # Set up logging
    logger = logging.getLogger('pycit_tester.py: %s' % time.asctime())
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
########################################################################
    # trimmomatic testing.
    logger.info("starting trimmomatic testing")
    trim = trimmomatic.Trimmomatic("trimmomatic")

    logger.info("Trim reads by quality")
    parameters = trimmomatic.Parameters(threads=4)
    steps = trimmomatic.Steps(ILLUMINACLIP="{0}:2:30:10".format(ADAPTERS))
    results = trim.run(LEFT_READS, RIGHT_READS, "trimmed_reads", PREFIX,
                       PHREDSCORE, parameters, steps)
    logger.info("Trimming returned:", results)
    logger.info("Trimming command: %s" % results.command)
    # get these exact file names from the named tuple
    logger.info("Trimming R1 file: %s" % results.outfileR1paired)
    logger.info("Trimming R2 file: %s" % results.outfileR2paired)
    logger.info("Trimming output: %s" % results.stderr)
    LEFT_TRIMMED = results.outfileR1paired
    RIGHT_TRIMMED = results.outfileR2paired

########################################################################
    # error correction testing
    error_correct = error_correction.Error_Correction("spades.py")
    logger.info("error correction using Bayes hammer")
    EC_results = error_correct.run(LEFT_TRIMMED, RIGHT_TRIMMED, THREADS,
                      "error_correction")
    # get these exact file names from the named tuple
    logger.info("ECorr R1 file: %s" % EC_results.Left_read_correct)
    logger.info("ECorr R2 file: %s" % EC_results.right_read_correct)
    logger.info("error correc output: %s" % EC_results.stderr)
    LEFT_TRIM_EC = EC_results.Left_read_correct
    RIGHT_TRIM_EC = EC_results.right_read_correct

########################################################################
    # PEAR testing - assemble
    assemble = pear.Pear("pear")
    logger.info("assembly using PEAR")
    results_pear = assemble.run(LEFT_TRIMMED, RIGHT_TRIMMED,
                 THREADS, "PEAR_assembled", PREFIX)

    logger.info("\n\nPEAR returned:", results_pear)
    logger.info("PEAR command: %s" % results_pear.command)
    logger.info("PEAR output: %s" % results_pear.stderr)
    PEAR_ASSEMBLED = results_pear.outfileassembled

    # Pear testing with EC reads
    logger.info("assembly using PEAR with error corrected reads")
    results_pear_EC = assemble.run(LEFT_TRIM_EC, RIGHT_TRIM_EC,
                 THREADS, "PEAR_assembled_EC", PREFIX)

    logger.info("PEAR returned using error corrected reads:", results_pear_EC)
    logger.info("PEAR command: %s" % results_pear_EC.command)
    logger.info("PEAR output: %s" % results_pear_EC.stderr)
    PEAR_ASSEMBLED_EC = results_pear_EC.outfileassembled

########################################################################
    # FLASH testing - assemble
    assemble = flash.Flash("flash")
    logger.info("assembly using Flash")
    results_flash = assemble.run(LEFT_TRIMMED, RIGHT_TRIMMED,
                                 THREADS, "flash_assembled", PREFIX)

    logger.info("\n\nFlash returned:", results_flash)
    logger.info("Flash command: %s" % results_flash.command)
    logger.info("Flash output: %s" % results_flash.stderr)
    FLASH_ASSEMBLED = results_flash.outfileextended

    # Flash testing with EC reads
    logger.info("assembly using Flash with error corrected reads")
    results_flash_EC = assemble.run(LEFT_TRIM_EC, RIGHT_TRIM_EC,
                                    THREADS, "flash_assembled_EC", PREFIX)

    logger.info("Flash returned using error corrected reads:", results_flash_EC)
    logger.info("Flash command: %s" % results_flash_EC.command)
    logger.info("Flash output: %s" % results_flash_EC.stderr)
    FLASH_ASSEMBLED_EC = results_flash_EC.outfileextended

########################################################################
    # cleaning up unwanted files
    logger.info("testing the removal of unwanted file")
    unwanted_list = ["unpaired_R*.fq.gz",
                     "PEAR_assembled/*.discarded*fastq",
                     "PEAR_assembled/*.unassembled*fastq",
                     "PEAR_assembled_EC/*.discarded*fastq",
                     "PEAR_assembled_EC/*.unassembled*fastq"]

    delete_file = clean_up.Clean_up("rm")
    for i in unwanted_list:
        try:
                del_results = delete_file.run(i)
                logger.info("removing file: %s" % i)
                logger.info("removing info: %s" % del_results)
        except:
                ValueError
                continue
        

########################################################################
    # convert format using seq_crumbs
    format_change = seq_crumbs.Convert_Format("convert_format",
                                              logger)
    format_change.run(PEAR_ASSEMBLED, "PEAR_assembled",
                      logger)
    format_change.run(PEAR_ASSEMBLED_EC, "PEAR_assembled_EC",
                      logger)
    format_change.run(FLASH_ASSEMBLED, "flash_assembled",
                      logger)
    format_change.run(FLASH_ASSEMBLED_EC, "flash_assembled_EC",
                      logger)

########################################################################
    # convert format using biopython
    
########################################################################
    # deduplicate read:
    #python deduplicate_rename.py
    #-f DNAMIX_S95_L001.assembled.fasta
    #-d database.out
    #-o temp.fasta
    logger.info("deduplicating reads")
    dedup_prog = os.path.join("pycits", "deduplicate_rename.py")
    dedup = deduplicate.Deduplicate(dedup_prog, logger)
    fasta_file = os.path.join("fasta_converted",
                              "DNAMIX_S95_L001.assembled.fasta")
    database = os.path.join("fasta_converted", "database.out")
    out_dedup = os.path.join("fasta_converted", "temp.fasta")
    dedup.run(dedup_prog, fasta_file, database, out_dedup,
              logger)

########################################################################
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
########################################################################
    # SWARM testing - assemble
    cluster = swarm.Swarm("swarm")
    assembled_fa_reads = os.path.join("fasta_converted",
                                      "temp.fasta")

    logger.info("clustering with Swarm")
    # def run(self, infname, outdir, threads, threshold, dry_run=False):
    cluster.run(assembled_fa_reads, "Swarm_cluster_assembled_reads_only",
                THREADS, 1)
    cluster_outdata = cluster.run(fasta_out, "Swarm_cluster", THREADS, 1, )

    logger.info("\n\nswarm returned using error corrected reads:", cluster_outdata)
    logger.info("swarm command: %s" % cluster_outdata.command)
    logger.info("swarm output: %s" % cluster_outdata.stderr)

########################################################################
    # recode cluster output
    #python pycits/parse_clusters_new_to_old_name.py -i
    # -i swarm_tempd1.out -d database.out -o decoded_clusters.out
    logger.info("renaming swarm output")
    rename_prog = os.path.join("pycits",
                               "parse_clusters_new_to_old_name.py")
    rename = rename_clusters_new_to_old.Rename(rename_prog, logger)
    infile = os.path.join("Swarm_cluster", "swarm.out")
    outfile = os.path.join("Swarm_cluster", 
               "swarm_clustering_d1_renamed.out")
    rename.run(rename_prog, infile, ITS_database,
               database, outfile, logger)

    # need to summarise the clusters now
    # how? scripts already made, are they good enough?

    #####################################################################
    post_analy_cmd = ["python",
           "./post_analysis/get_results_from_cluster_and_novel_clusterings.py",
            "-f", "fasta_converted/DNAMIX_S95_L001.assembled.fasta",
            "--all_fasta", "fasta_converted/database_assembled_reads.fasta",
            "--seq_db", "data/ITS_database_NOT_confirmed_correct_last14bases_removed.fasta",
            "--min_novel_cluster_threshold", "2",
            "--left", "./data/reads/DNAMIX_S95_L001_R1_001.fastq",
            "--right", "./data/reads/DNAMIX_S95_L001_R2_001.fastq",
            "--Name_of_project",  "testing_output",
            "--in",  "Swarm_cluster/swarm.out",
            "--difference", "1",
            "-o", "testing_script",
            "--old_to_new", "fasta_converted/database.out"]
    post_analy_cmd = ' '.join(post_analy_cmd)
    print (post_analy_cmd)
    pipe = subprocess.run(post_analy_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    
