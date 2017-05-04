#!/usr/bin/env python
# THIS IS currently under development.
#
# TODO - put this functionality under a proper test framework (e.g. nosetests
#        or equivalent) - in progress
# All modules used are subjected to nosetests

import sys
import errno
import logging
import logging.handlers
import multiprocessing
import os
import time
import traceback
import subprocess
import sklearn
import argparse
import pysam

from pycits.tools import convert_fq_to_fa, NotExecutableError, trim_seq,\
     dereplicate_name, check_OTU_db_abundance_val, \
     parse_tab_file_get_clusters, filter_sam_file, reformat_cdhit_clustrs,\
     reformat_sam_clusters, reformat_swarm_cls, reformat_blast6_clusters

from pycits.Rand_index import pairwise_comparison_Rand

from pycits import tools, fastqc, trimmomatic, pear, error_correction,\
    flash, clean_up, swarm, seq_crumbs, bowtie_build, bowtie_map,\
    cd_hit, blast, vsearch, samtools_index, muscle

# check this is running under python 3
if sys.version_info[0] < 3:
    # break the program
    raise ImportError ("\n\t\tYou must be using Python 3\n")
    sys.exit(1)

if "-v" in sys.argv or "--version" in sys.argv:
    print("Pycits classify OTU: v0.0.1")
    sys.exit()

#########################################################################


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding " +
                                     "This is currently a draft.  " ,
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("METAPY")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--thread", dest='threads',
                          action="store", default="4",
                          type=str,
                          help="number of threads")

    optional.add_argument("-l", "--left", dest='left',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "reads",
                                               "DNAMIX_S95_" +
                                               "L001_R1_001.fastq"),
                          type=str,
                          help="left illumina reads, " +
                          "default is for tests")

    optional.add_argument("-r", "--right", dest='right',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "reads",
                                               "DNAMIX_S95_" +
                                               "L001_R2_001.fastq"),
                          type=str,
                          help="right illumina reads")

    optional.add_argument("--OTU_DB", dest='OTU_DB',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "ITS_db_NOT_conf_correct" +
                                               "_last14bp_removd.fasta"),
                          type=str,
                          help="right illumina reads")

    optional.add_argument("-a", "--assemble", dest='assemble',
                          action="store",
                          choices=["pear", "flash"],
                          help="program to assemble with " +
                          "flash or pear " +
                          "default: %(default)s", default="pear",
                          type=str)

    optional.add_argument("--adaptors", dest='adaptors',
                          action="store",
                          default= os.path.join(file_directory,
                                                "adapters",
                                                "TruSeq3-PE.fa"),
                          type=str,
                          help="adaptors for trimming. Can supply custom " +
                          " file if desired")

    optional.add_argument("--left_trim", dest='left_trim',
                          action="store",
                          default=53,
                          type=int,
                          help="left_trim for primers or conserved " +
                          "regions. Default 53 ")

    optional.add_argument("--right_trim", dest='right_trim',
                          action="store",
                          default= 0,
                          type=int,
                          help="right_trim for primers or conserved " +
                          "regions.  Default 0")

    optional.add_argument("--phred", dest='phred',
                          action="store",
                          default="phred33",
                          type=str,
                          help="phred33 is default. " +
                          "Dont change unless sure")

    optional.add_argument("--cdhit_threshold", dest='cdhit_threshold',
                          action="store",
                          default="0.99",
                          type=str,
                          help="percentage identify for cd-hit " +
                          "Default -0.99")

    optional.add_argument("--swarm_d_value", dest='swarm_d_value',
                          action="store",
                          default=1,
                          type=int,
                          help="the difference d value for clustering " +
                          "in swarm. Default 1")

    optional.add_argument("--blastclust_threshold", dest='blastclust_threshold',
                          action="store",
                          default=0.90,
                          type=float,
                          help="the threshold for blastclust clustering " +
                          " Default -S 0.90")

    optional.add_argument("--vesearch_threshold", dest='vesearch_threshold',
                          action="store",
                          default=0.99,
                          type=float,
                          help="the threshold for vsearch clustering " +
                          " Default 0.99")

    optional.add_argument("--verbose", dest="verbose",
                          action="store_true",
                          default=False,
                          help="Report verbose output")

    optional.add_argument("--Error_correction",
                          dest="Error_correction",
                          action="store_true",
                          default=True,
                          help="to perform Illumina error correction ")

    optional.add_argument("--align", dest="align",
                          action="store_true",
                          default=False,
                          help="to align clusters in the output " +
                          "you must have muscle in your PATH as muscle")

    optional.add_argument("--percent_identity", dest="percent_identity",
                          action="store_true",
                          default=False,
                          help="blast the cluster to return " +
                          "pairwise percentage identity")

    optional.add_argument("--min_novel_cluster_threshold",
                          dest="min_novel_cluster_threshold",
                          type=str,
                          default="2",
                          help="min size of a cluster to consider as real " +
                          "anything smaller than this is ignored")

    optional.add_argument("--logfile", dest="logfile",
                          action="store", default="pipeline.log",
                          type=str,
                          help="Logfile name")

    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message"
                          " type --version for version")
    args = parser.parse_args()
    return args, file_directory


def decompress(infile):
    """function to decompress gzipped reads"""
    cmd = ["gunzip",
           infile]
    cmd = ' '.join(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    return infile.split(".gz")[0]


def compress(infile):
    """function to decompress reads, make them .gz"""
    cmd = ["gzip",
           infile]
    cmd = ' '.join(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

###################################################################
# Global variables
WARNINGS = ""
args, FILE_DIRECTORY = get_args()
# setting up some test variables
THREADS = args.threads
# test read set are defaults
LEFT_READS = args.left
if LEFT_READS.endswith(".gz"):
    LEFT_READS = decompress(LEFT_READS)
RIGHT_READS = args.right
if RIGHT_READS.endswith(".gz"):
    RIGHT_READS = decompress(RIGHT_READS)
assert(LEFT_READS != RIGHT_READS)
READ_PREFIX = os.path.split(LEFT_READS)[-1].split("_R")[0]
PHREDSCORE = args.phred
PREFIX = READ_PREFIX
ADAPTERS = args.adaptors
OUTDIR_TRIM = "trimmed_reads"
OUTFILES = [os.path.join(OUTDIR_TRIM, PREFIX + suffix) for suffix in
           ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
            "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]
OTU_DATABASE = args.OTU_DB
WORKING_DIR = os.getcwd()
# set this to false for now to not run it
SEQ_CRUMBS = False
CDHIT_THRESHOLD = str(args.cdhit_threshold)
if float(CDHIT_THRESHOLD) < 0.8:
    sys.exit("\nCDHIT_THRESHOLD must be less than 0.8\n")
SWARM_D_VALUE = args.swarm_d_value
if SWARM_D_VALUE > 2:
    WARNINGS = "swarm threshold is very 'loose' = %d\n" % SWARM_D_VALUE
VSEARCH_THRESHOLD = args.vesearch_threshold
ERROR_CORRECTION = args.Error_correction
ASSEMBLE_PROG = args.assemble

# left_trim 53 - this should be default??
LEFT_TRIM = args.left_trim
RIGHT_TRIM = args.right_trim
assert(READ_PREFIX == os.path.split(RIGHT_READS)[-1].split("_R")[0])

CLUSTER_FILES_FOR_RAND_INDEX = []
RESULTS = []
#####################################################################


def check_tools_exist():
    """function to check to see what tools are in the PATH,
    decide what we can use
    Returns a list of programs that were exectable and a warning string.
    The warnings are tools that were not executable"""
    tools_list = []
    Warning_out = WARNINGS + "Tool executable warning: "
    try:
        flash.Flash("flash")
        tools_list.append("flash")
    except ValueError:
        Warning_out = Warning_out + "Flash not in path"
    try:
        error_correction.Error_Correction("spades.py")
        tools_list.append("error_correction")
    except ValueError:
        Warning_out = Warning_out + "spades.py not in path\n"
    try:
        vsearch.Vsearch_derep("vsearch")
        tools_list.append("vsearch")
    except ValueError:
        Warning_out = Warning_out + "vsearch not in path\n"
    try:
        trimmomatic.Trimmomatic("trimmomatic")
        tools_list.append("trimmomatic")
    except ValueError:
        Warning_out = Warning_out + "trimmomatic not in path\n"
    try:
        swarm.Swarm("swarm")
        tools_list.append("swarm")
    except ValueError:
        Warning_out = Warning_out + "swarm not in path\n"
    try:
        samtools_index.Samtools_Index("samtools")
        tools_list.append("samtools")
    except ValueError:
        Warning_out = Warning_out + "samtools not in path\n"
    try:
        pear.Pear("pear")
        tools_list.append("pear")
    except ValueError:
        Warning_out = Warning_out + "pear not in path\n"
    try:
        muscle.Muscle("muscle")
        tools_list.append("muscle")
    except ValueError:
        Warning_out = Warning_out + "muscle not in path\n"
    try:
        fastqc.FastQC("fastqc")
        tools_list.append("fastqc")
    except ValueError:
        Warning_out = Warning_out + "fastqc not in path\n"
    try:
        fastqc.FastQC("cd-hit-est")
        tools_list.append("cd-hit-est")
    except ValueError:
        Warning_out = Warning_out + "cd-hit-est not in path\n"
    try:
        bowtie_build.Bowtie2_Build("bowtie2-build")
        tools_list.append("bowtie2")
    except ValueError:
        Warning_out = Warning_out + "bowtie2 not in path\n"
    try:
        blast.Blastclust("blastclust")
        tools_list.append("blastclust")
    except ValueError:
        Warning_out = Warning_out + "blastclust not in path\n"
    return tools_list, Warning_out


def make_folder(folder):
    """function to make a folder with desired name"""
    dest_dir = os.path.join(WORKING_DIR, folder)
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists " +
               "I will write over what is in there!!")
    return dest_dir


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type,
                                              exc_value,
                                              exc_traceback))


def covert_chop_read(infile):
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
    logger.info("chopping the reads %s" % infile)
    trim_seq(infile + ".bio.fasta",
             infile + ".bio.chopped.fasta",
             LEFT_TRIM, RIGHT_TRIM)

#######################################################################
# Run as script
if __name__ == '__main__':
    # Set up logging
    logger = logging.getLogger('METAPY.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(args.logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        outstr = "Could not open %s for logging" % args.logfile
        logger.error(outstr)
        sys.exit(1)
    # Report input arguments
    outstr = "Command-line: %s" % ' '.join(sys.argv)
    logger.info(outstr)
    outstr = "Starting testing: %s" % time.asctime()
    logger.info(outstr)
    # Get a list of tools in path!
    logger.info("checking which programs are in PATH")
    tools_list, Warning_out = check_tools_exist()
    logger.info(Warning_out)

    ####################################################################
    # fastqc QC of raw reads
    if "fastqc" in tools_list:
        FASTQC_FOLDER = make_folder(PREFIX + "_fastqc")
        logger.info("starting fastqc")
        outstr = "made folder  %s" % FASTQC_FOLDER
        logger.info(outstr)
        qc = fastqc.FastQC("fastqc")
        qc_results = qc.run(LEFT_READS, PREFIX + FASTQC_FOLDER)
        outstr = "fastqc output: %s" % qc_results.command
        logger.info(outstr)
        outstr = "fastqc stderr: %s" % qc_results.stderr
        logger.info(outstr)

    ####################################################################
    # trimmomatic trm reads
    if "trimmomatic" in tools_list:
        TRIM_FOLDER = make_folder(PREFIX + "_timmomatic")

        logger.info("starting trimmomatic testing")
        outstr = "made folder %s" % TRIM_FOLDER
        logger.info(outstr)
        trim = trimmomatic.Trimmomatic("trimmomatic")

        logger.info("Trim reads by quality")
        parameters = trimmomatic.Parameters(threads=4)
        steps = trimmomatic.Steps(ILLUMINACLIP="{0}:2:30:10".format(ADAPTERS))
        results = trim.run(LEFT_READS, RIGHT_READS,
                           TRIM_FOLDER, PREFIX,
                           PHREDSCORE, parameters,
                           steps)
        logger.info("Trimming returned:", results)
        outstr = "Trimming command: %s" % results.command
        logger.info(outstr)
        # get these exact file names from the named tuple
        outstr = "Trimming output: %s" % results.stderr
        logger.info(outstr)
        LEFT_TRIMMED = results.outfileR1paired
        RIGHT_TRIMMED = results.outfileR2paired

    ####################################################################
    # error correction testing
    if "error_correction" in tools_list:
        if ERROR_CORRECTION == True:
            # we will error correct the read and reasign LEFT_TRIMMED
            # with the EC reads
            EC_FOLDER = make_folder(PREFIX + "_EC")
            error_corr = error_correction.Error_Correction("spades.py")
            logger.info("error correction using Bayes hammer")
            outstr = "made folder %s" % EC_FOLDER
            logger.info(outstr)
            EC_results = error_corr.run(LEFT_TRIMMED,
                                        RIGHT_TRIMMED,
                                        THREADS,
                                        EC_FOLDER)
            # get these exact file names from the named tuple
            outstr = "error correc output: %s" % EC_results.stderr
            logger.info(outstr)
            # assign the EC to the LEFT and RIGHT trimmed variable
            LEFT_TRIMMED = EC_results.Left_read_correct
            RIGHT_TRIMMED = EC_results.right_read_correct
            logger.info("Trimmed reads have been error corrected " +
                        "using Bayes hammer")

    ####################################################################
    # PEAR testing - assemble
    if "pear" in tools_list and ASSEMBLE_PROG == "pear":
        logger.info("assembly using PEAR")
        if ERROR_CORRECTION == True:
            suffix = "_PEAR_EC"
            logger.info("PEAR will use trimmed and error corrected reads")
        else:
            suffix = "_PEAR"
        ASSEMBLY_FOLDER = make_folder(PREFIX + suffix)
        # call the class
        assemble = pear.Pear("pear")
        results_pear = assemble.run(LEFT_TRIMMED,
                                    RIGHT_TRIMMED,
                                    THREADS,
                                    ASSEMBLY_FOLDER,
                                    PREFIX)
        logger.info("\n\nPEAR returned:", results_pear)
        outstr = "PEAR command: %s" % results_pear.command
        logger.info(outstr)
        outstr = "PEAR output: %s" % results_pear.stderr
        logger.info(outstr)
        ASSEMBLED = results_pear.outfileassembled
        # call the function
        covert_chop_read(ASSEMBLED)

    ####################################################################
    # FLASH testing - assemble
    if "flash" in tools_list and ASSEMBLE_PROG == "flash":
        logger.info("assembly using FLASH")
        if ERROR_CORRECTION == True:
            suffix = "_FLASH_EC"
            logger.info("FLASH will use trimmed and error corrected reads")
        else:
            suffix = "_FLASH"
        ASSEMBLY_FOLDER = make_folder(PREFIX + suffix)
        assemble = flash.Flash("flash")
        logger.info("assembly using Flash")
        results_flash = assemble.run(LEFT_TRIMMED,
                                     RIGHT_TRIMMED,
                                     THREADS,
                                     ASSEMBLY_FOLDER,
                                     PREFIX)
        logger.info("\n\nFlash returned:", results_flash)
        outstr = "Flash command: %s" % results_flash.command
        logger.info(outstr)
        outstr = "Flash output: %s" % results_flash.stderr
        logger.info(outstr)
        ASSEMBLED = results_flash.outfileextended
        # call the function
        covert_chop_read(ASSEMBLED)

    ####################################################################
    # convert format using seq_crumbs
    if SEQ_CRUMBS == True:
        format_change = seq_crumbs.Convert_Format("convert_format",
                                                  logger)
        format_change.run(ASSEMBLED,
                          ASSEMBLY_FOLDER,
                          logger)
    # first cat the db and EC, trimmed reads.
    cat_cmd = ["cat", OTU_DATABASE,
               ASSEMBLED + ".bio.chopped.fasta",
               ">",
               "assembled_fa_and_OTU_db.fasta"]
    cat_cmd = ' '.join(cat_cmd)
    outstr = "I going to cat these files %s" % cat_cmd
    logger.info(outstr)
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

    ####################################################################
    # deduplicate reads:
    # use the function in tools dereplicate_name()
    if "swarm" in tools_list:
        logger.info("deduplicating reads")
        # dereplicate_name (infasta, db_out, out_fasta)
        dereplicate_name(ASSEMBLED + ".bio.chopped.fasta",
                         "db_old_to_new_names.txt",
                         ASSEMBLED + "for_swarm.fasta")

    ####################################################################
    # SWARM testing - assemble
    if "swarm" in tools_list:
        swarm_parameters = swarm.Parameters(t=1, d=SWARM_D_VALUE)
        SWARM_FOLDER = make_folder(PREFIX + "_Swarm_d%d" % SWARM_D_VALUE)
        cluster = swarm.Swarm("swarm")
        assembled_fa_reads = ASSEMBLED + "for_swarm.fasta"
        logger.info("clustering with Swarm")
        # need to check the OTU database has abundance value
        # call the function from tools
        outstr = "OTU was %s" % OTU_DATABASE
        logger.info(outstr)
        OTU_DATABASE_SWARM = check_OTU_db_abundance_val(OTU_DATABASE)
        outstr = "OTU is %s" % OTU_DATABASE_SWARM
        logger.info(outstr)
        # need to cat the assembled_fasta with the database
        cat_cmd = ["cat",
                   OTU_DATABASE_SWARM,
                   assembled_fa_reads,
                   ">",
                   "assembled_reads_and_OTU_db.fasta"]
        cat_cmd = ' '.join(cat_cmd)
        outstr = "I cat these files %s" % cat_cmd
        logger.info(outstr)
        pipe = subprocess.run(cat_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        cluster_outdata = cluster.run("assembled_reads_and_OTU_db.fasta",
                                      SWARM_FOLDER,
                                      swarm_parameters)
        # use named tuple to get the outpfile name
        SWARM_OUT = cluster_outdata.outfilename

        logger.info("swarm returned using error corrected reads:",
                    cluster_outdata)
        outstr = "swarm command: %s" % cluster_outdata.command
        logger.info(outstr)
        outstr = "swarm output: %s" % cluster_outdata.stderr
        logger.info(outstr)

        ################################################################
        # recode cluster output
        logger.info("renaming swarm output")
        # parse_tab_file_get_clusters(filename1, database, out_file)
        parse_tab_file_get_clusters(SWARM_OUT,
                                    OTU_DATABASE_SWARM,
                                    "db_old_to_new_names.txt",
                                    SWARM_OUT + "RENAMED_abundance",
                                    True)
        parse_tab_file_get_clusters(SWARM_OUT,
                                    OTU_DATABASE,
                                    "db_old_to_new_names.txt",
                                    SWARM_OUT + "RENAMED")
        # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
        logger.info("reformatting swarm output for post analysis")

        reformat_swarm_cls(SWARM_OUT + "RENAMED",
                           OTU_DATABASE,
                           "assembled_fa_and_OTU_db.fasta",
                           SWARM_OUT + "for_R",
                           False)
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(SWARM_OUT + "for_R")
        cmd_s = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_reads_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE_SWARM,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(SWARM_FOLDER, "clusters"),
                 "--in",
                 SWARM_OUT + "RENAMED_abundance",
                 "--difference", str(SWARM_D_VALUE),
                 "-o",
                 os.path.join(SWARM_FOLDER,
                              "%s_swarm_results_%d.RESULTS" %
                              (PREFIX, SWARM_D_VALUE)),
                 "--old_to_new", "db_old_to_new_names.txt"]
        swarm_result = os.path.join(SWARM_FOLDER,
                                    "%s_swarm_results_%d.RESULTS" %
                                    (PREFIX, SWARM_D_VALUE))
        RESULTS.append("swarm\t%s" % swarm_result)
        cmd_s = ' '.join(cmd_s)
        if args.align == True:
            logger.info("going to align the cluster. Will take ages!")
            cmd_s = cmd_s + " --align True"
        if args.percent_identity == True:
            cmd_s = cmd_s + " --blast True"
        outstr = "%s = post analysis command" % cmd_s
        logger.info(outstr)
        pipe = subprocess.run(cmd_s, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        logger.info("graphically represent swarm clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    SWARM_OUT + "RENAMED_abundance"
                    " --db",
                    OTU_DATABASE_SWARM]
        plot_cmd = ' '.join(plot_cmd)
        logger.info("plotting command = %s" % plot_cmd)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    #####################################################################
    # run cd hit
    # first cat the db and EC, trimmed reads.
    cat_cmd = ["cat", OTU_DATABASE,
               ASSEMBLED + ".bio.chopped.fasta",
               ">",
               "assembled_fa_and_OTU_db.fasta"]
    cat_cmd = ' '.join(cat_cmd)
    outstr = "I cat these files %s" % cat_cmd
    logger.info(outstr)
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    if "cd-hit-est" in tools_list:
        CDHIT_FOLDER = make_folder(PREFIX + "_cd_hit_%s" % CDHIT_THRESHOLD)
        cluster = cd_hit.Cd_hit("cd-hit-est")
        results = cluster.run("assembled_fa_and_OTU_db.fasta",
                              THREADS,
                              CDHIT_THRESHOLD,
                              CDHIT_FOLDER,
                              PREFIX)
        outstr = "cdhit: %s" % results.command
        logger.info(outstr)
        logger.info("reformatting cd hit output")
        reformat_cdhit_clustrs(results.clusters,
                               results.clusters + "_1_line_per",
                               results.clusters + "for_R")
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(results.clusters + "for_R")
        # analyse the clusters
        cmd_c = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(CDHIT_FOLDER, "clusters"),
                 "--in",
                 results.clusters + "_1_line_per",
                 "--difference", CDHIT_THRESHOLD,
                 "-o",
                 os.path.join(CDHIT_FOLDER,
                              "%s_cdhit_results_%s.RESULTS" %
                              (PREFIX, str(CDHIT_THRESHOLD))),
                 "--old_to_new", "db_old_to_new_names.txt"]

        cd_hit_result = os.path.join(CDHIT_FOLDER,
                                     "%s_cdhit_results_%s.RESULTS" %
                                     (PREFIX, str(CDHIT_THRESHOLD)))
        RESULTS.append("cdhit\t%s" % cd_hit_result)

        cmd_c = ' '.join(cmd_c)
        if args.align == True:
            logger.info("going to align the cluster. Will take ages!")
            cmd_c = cmd_c + " --align True"
        if args.percent_identity == True:
            cmd_c = cmd_c + " --blast True"
        outstr = "%s = post analysis command" % cmd_c
        logger.info(outstr)
        pipe = subprocess.run(cmd_c, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        logger.info("graphically represent cdhit clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    results.clusters + "_1_line_per",
                    " --db",
                    OTU_DATABASE]
        plot_cmd = ' '.join(plot_cmd)
        logger.info("plotting command = %s" % plot_cmd)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)


    ####################################################################
    # run vsearch
    if "vsearch" in tools_list:
        VSEARCH_FOLDER = make_folder(PREFIX + "_Vsearch_%.2f" %
                                     VSEARCH_THRESHOLD)

        V_derep = vsearch.Vsearch_derep("vsearch")
        derep_results = V_derep.run(ASSEMBLED + ".bio.chopped.fasta",
                                    VSEARCH_FOLDER,
                                    PREFIX)
        outstr = "vsearch derep: %s" % derep_results.command
        logger.info(outstr)
        logger.info("vsearch clustering")
        V_clst = vsearch.Vsearch_cluster("vsearch")
        vclu_rlts = V_clst.run(derep_results.fasta,
                               VSEARCH_FOLDER,
                               PREFIX,
                               OTU_DATABASE,
                               THREADS,
                               VSEARCH_THRESHOLD)
        outstr = "vsearch cluster: %s" % vclu_rlts.command
        logger.info(outstr)
        reformat_blast6_clusters(vclu_rlts.blast6,
                                 "assembled_fa_and_OTU_db.fasta",
                                 vclu_rlts.blast6 + "for_R")
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(vclu_rlts.blast6 + "for_R")

        cmd_v = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(VSEARCH_FOLDER, "clusters"),
                 "--in",
                 vclu_rlts.blast6 + "for_R_1_line",
                 "--difference", str(VSEARCH_THRESHOLD),
                 "-o",
                 os.path.join(VSEARCH_FOLDER,
                              "%s_vsearch_results_%s.RESULTS" %
                              (PREFIX, str(VSEARCH_THRESHOLD))),
                 "--old_to_new", "db_old_to_new_names.txt"]

        v_result = os.path.join(VSEARCH_FOLDER,
                                "%s_vsearch_results_%s.RESULTS" %
                                (PREFIX, str(VSEARCH_THRESHOLD)))
        RESULTS.append("vsearch\t%s" % v_result)

        cmd_v = ' '.join(cmd_v)
        if args.align == True:
            logger.info("going to align the cluster. Will take ages!")
            cmd_v = cmd_v + " --align True"
        if args.percent_identity == True:
            cmd_v = cmd_v + " --blast True"
        outstr = "%s = post analysis command" % cmd_v
        logger.info(outstr)
        pipe = subprocess.run(cmd_v, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        logger.info("graphically represent vsearch default clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    vclu_rlts.blast6 + "for_R_1_line",
                    " --db",
                    OTU_DATABASE]
        plot_cmd = ' '.join(plot_cmd)
        outstr = "plotting command = %s" % plot_cmd
        logger.info(outstr)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    ###################################################################
    # deduplicate reads:
    # use the function in tools dereplicate_name()
    # run vsearch
    if "vsearch" in tools_list:
        logger.info("deduplicating reads")
        # dereplicate_name (infasta, db_out, out_fasta)
        # derep the databse for vsearch
        derep_db = V_derep.run(OTU_DATABASE,
                               VSEARCH_FOLDER,
                               PREFIX)
        dereplicate_name(ASSEMBLED + ".bio.chopped.fasta",
                         "db_old_to_new_names_vsearch.txt",
                         ASSEMBLED + "for_vsearch.fasta",
                         True)
        # cat the derep reads and derep db together

        cat_cmd = ["cat", derep_db.fasta,
                   ASSEMBLED + "for_vsearch.fasta",
                   ">",
                   "assembled_fa_and_OTU_db_vesearch.fasta"]
        cat_cmd = ' '.join(cat_cmd)
        outstr = "I cat these files %s" % cat_cmd
        logger.info(outstr)
        pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

        # use vsearch to get concensus and aligned clusters
        clu_fasta = vsearch.Vsearch_fastas("vsearch")
        fasta_results = clu_fasta.run("assembled_fa_and_OTU_db" +
                                      "_vesearch.fasta",
                                      VSEARCH_FOLDER,
                                      PREFIX,
                                      OTU_DATABASE,
                                      THREADS,
                                      VSEARCH_THRESHOLD)

        reformat_blast6_clusters(fasta_results.blast6,
                                 "assembled_fa_and_OTU_db.fasta",
                                 fasta_results.blast6 + "for_R")
        outstr = "vsearch cluster: %s" % vclu_rlts.command
        logger.info(outstr)
        CLUSTER_FILES_FOR_RAND_INDEX.append(fasta_results.blast6 + "for_R")

        cmd_F = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(VSEARCH_FOLDER, "clusters_clusterfast"),
                 "--in",
                 fasta_results.blast6 + "for_R_1_line",
                 "--difference", str(VSEARCH_THRESHOLD),
                 "-o",
                 os.path.join(VSEARCH_FOLDER,
                              "%s_vserach_clu_fasta_results_%s.RESULTS" %
                              (PREFIX, str(VSEARCH_THRESHOLD))),
                 "--old_to_new",
                 "db_old_to_new_names.txt"]

        v_f_result = os.path.join(VSEARCH_FOLDER,
                                  "%s_vserach_clu_fasta_results_%s.RESULTS" %
                                  (PREFIX, str(VSEARCH_THRESHOLD)))

        RESULTS.append("vse_faclus\t%s" % v_f_result)

        cmd_F = ' '.join(cmd_F)
        if args.align == True:
            logger.info("going to align the cluster. Will take ages!")
            cmd_F = cmd_F + " --align True"
        if args.percent_identity == True:
            cmd_F = cmd_F + " --blast True"
        outstr = "%s = post analysis command" % cmd_v
        logger.info(outstr)
        pipe = subprocess.run(cmd_F, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        logger.info("graphically represent vsearch clu_fasta clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    fasta_results.blast6 + "for_R_1_line",
                    " --db",
                    OTU_DATABASE]
        plot_cmd = ' '.join(plot_cmd)
        outstr = "plotting command = %s" % plot_cmd
        logger.info(outstr)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    #####################################################################
    # MAP THE READS WITH BOWTIE
    if "bowtie2" in tools_list:
        BOWTIE_FOLDER = make_folder(PREFIX + "_bowtie")
        obj = bowtie_build.Bowtie2_Build("bowtie2-build")
        results = obj.run(OTU_DATABASE, "OTU")
        outstr = "bowtie build %s" % results.command
        logger.info(outstr)
        obj = bowtie_map.Bowtie2_Map("bowtie2")
        results = obj.run(ASSEMBLED + ".bio.chopped.fasta",
                          "OTU",
                          BOWTIE_FOLDER,
                          THREADS)
        outstr = "bowtie map %s" % results.command
        logger.info(outstr)
        logger.info("pysam to filter the mapping")
        samfile = pysam.AlignmentFile(results.sam, "r")

        # call the function from tools
        # this return matches with zero mismatches, but not to be interpreted
        # as a perfect match?!?!
        cig_list, matches = filter_sam_file(results.sam,
                                            (os.path.join(BOWTIE_FOLDER,
                                                          "pysam_perfect_" +
                                                          "cigar.txt")))
        outstr = "pysam found %d perfect(?) cigar MATCHES" % len(cig_list)
        logger.info(outstr)
        # using grep to get perfect matches:
        grep_cmd = ' '.join(['cat',
                             results.sam,
                             '|',
                             'grep',
                             '"AS:i:0"',
                             '>',
                             results.sam + "perfect_map"])
        outstr = "grep for perfect reads %s" % grep_cmd
        logger.info(outstr)
        pipe = subprocess.run(grep_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        # sort these and make them unique.
        # for info only, this next command is not needed
        grep_cmd = ' '.join(['cat',
                             results.sam,
                             '|',
                             'grep',
                             '"AS:i:0"',
                             '|',
                             'cut -f3',
                             '|',
                             'sort'
                             '|',
                             'uniq'
                             '>',
                             results.sam + "perfect_map_name"])
        outstr = "grep for perfect reads %s" % grep_cmd
        logger.info(outstr)
        pipe = subprocess.run(grep_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        perfect_count = 0
        with open(results.sam + "perfect_map_name") as perfect_names:
            data = perfect_names.read().split("\n")
            for name in data:
                perfect_count = perfect_count + 1
                logger.info("%s = perfect match" % name)
        logger.info("%d = number perfect match" % perfect_count)
        reformat_sam_clusters(results.sam + "perfect_map",
                              "assembled_fa_and_OTU_db.fasta",
                              results.sam +
                              "for_R")
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(results.sam + "for_R")

        cmd_b = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(BOWTIE_FOLDER, "clusters_perfect_map"),
                 "--in",
                 results.sam + "for_R_1_line",
                 "--difference", "perfect_map",
                 "-o",
                 os.path.join(BOWTIE_FOLDER,
                              "bowtie_perfect.RESULTS"),
                 "--old_to_new", "db_old_to_new_names.txt"]
        bow_result = os.path.join(BOWTIE_FOLDER,
                                  "bowtie_perfect.RESULTS")

        RESULTS.append("bowtie\t%s" % bow_result)

        cmd_b = ' '.join(cmd_b)
        logger.info("%s = post analysis command" % cmd_b)
        pipe = subprocess.run(cmd_b, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        logger.info("graphically represent bowtie clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    results.sam + "for_R_1_line",
                    " --db",
                    OTU_DATABASE]
        plot_cmd = ' '.join(plot_cmd)
        outstr = "plotting command = %s" % plot_cmd
        logger.info(outstr)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    #####################################################################
    # run BLASTCLUST
    if "blastclust" in tools_list:
        blastclust_threshold = args.blastclust_threshold  # for now
        BLASTCL_FOLDER = make_folder(PREFIX + "_blastclust_%s" %
                                     str(blastclust_threshold))
        logger.info("running blastclust")
        bc = blast.Blastclust("blastclust")
        result = bc.run("assembled_fa_and_OTU_db.fasta", BLASTCL_FOLDER,
                        THREADS)

        # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
        # use this to reformat the blastclust clusters
        result_file_r = (os.path.join(BLASTCL_FOLDER,
                                      "assembled_fa_and_OTU_db.fasta" +
                                      ".blastclust99.forR"))
        reformat_swarm_cls(os.path.join(BLASTCL_FOLDER,
                                        "assembled_fa_and_OTU_db.fasta" +
                                        ".blastclust99.lst"),
                           OTU_DATABASE,
                           "assembled_fa_and_OTU_db.fasta",
                           result_file_r,
                           False)
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(result_file_r)
        logger.info("graphically represent swarm clusters")
        plot_cmd = ["python",
                    os.path.join(FILE_DIRECTORY,
                                 "bin",
                                 "draw_bar_chart_of_clusters.py"),
                    "-i",
                    result.outfilename,
                    " --db",
                    OTU_DATABASE]
        plot_cmd = ' '.join(plot_cmd)
        outstr = "plotting command = %s" % plot_cmd
        logger.info(outstr)
        pipe = subprocess.run(plot_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        cmd_b = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(BLASTCL_FOLDER, "clusters_%s" %
                              str(blastclust_threshold)),
                 "--in",
                 result.outfilename,
                 "--difference", str(blastclust_threshold),
                 "-o",
                 os.path.join(BLASTCL_FOLDER,
                              "%s_blastclust_results_%s.RESULTS" %
                              (PREFIX, str(blastclust_threshold))),
                 "--old_to_new", "db_old_to_new_names.txt"]

        bc_result = os.path.join(BLASTCL_FOLDER,
                                "%s_blastclust_results_%s.RESULTS" %
                                (PREFIX, str(blastclust_threshold)))
        RESULTS.append("blastclust\t%s" % bc_result)

        cmd_b = ' '.join(cmd_b)
        if args.align == True:
            logger.info("going to align the cluster. Will take ages!")
            cmd_b = cmd_b + " --align True"
        if args.percent_identity == True:
            cmd_b = cmd_b + " --blast True"
        outstr = "%s = post analysis command" % cmd_b
        logger.info(outstr)
        pipe = subprocess.run(cmd_b, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    Rand_results = pairwise_comparison_Rand(CLUSTER_FILES_FOR_RAND_INDEX,
                                            PREFIX +
                                            "_Rand_comparison.txt")
    for comp in result:
        outstr = "Rand comparison: %s" % comp
        logger.info(outstr)

    # compress the reads to save space
    # compress(LEFT_READS)
    # compress(RIGHT_READS)
    logger.info("compressed the reads")

    #####################################################################
    # compare the results files
    compare_prog = os.path.join(FILE_DIRECTORY,
                                "bin",
                                "compare_results.py")
    blastclust_threshold = str(blastclust_threshold)
    comp = "%s_RESULTS_cd_%s_sw_%s_BC_%s_V_%s.txt" % (PREFIX,
                                                      CDHIT_THRESHOLD,
                                                      str(SWARM_D_VALUE),
                                                      blastclust_threshold,
                                                      str(VSEARCH_THRESHOLD))
    # write the result file name to file. Easier to get these for the
    # next part - compare_prog
    f_out = open("temp.txt", "w")
    for i in RESULTS:
        f_out.write(i + "\n")
    f_out.close()
    cmd_r = " ".join(["python",
                      compare_prog,
                      " -o",
                      comp,
                      " --in_list",
                      "temp.txt"])
    outstr = "%s = comparison comment" % cmd_r
    logger.info(outstr)
    pipe = subprocess.run(cmd_r, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    outstr = "Pipeline complete: %s" % time.asctime()
    logger.info(outstr)
    os.remove("temp.txt")
