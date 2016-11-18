#!/usr/bin/env python
#
# Tools for working with trimmomatic
#
# trimmomatic:
# http://www.bioinformatics.babraham.ac.uk/projects/trimmomatic/
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe


import os
import sys

import subprocess
from .tools import is_exe


class Trimmomatic(object):
    """Class for working with trimmomatic"""
    def __init__(self, exe_path, logger=False):
        """Instantiate with location of executable"""
        if logger:
            self._logger = logger
        self._no_run = False
        if not os.path.isfile(exe_path):
            # TO DO: return this error message when it fails
            msg = ''.join(["trimmomatic is not valid trimming of "
                            "reads will not be run...           \n",
                            " SOLUTION: put trimmomatic in PATH   ",
                            " and this is still failing you may   ",
                            "need to rename/ copy your trimmomatic ",
                            "binary to a file called timmomatic. ",
                            " please make ensure file is executable ",
                            " and in your PATH  you can download ",
                            " the binaries from http://www.usadellab",
                            ".org/cms/?page=trimmomatic"])
            if logger:
                self._logger.warning(msg)
            self._no_run = True
        self._exe_path = exe_path

    def __build_cmd(self, trimmo_prog, L_reads, R_reads, threads,
                    outdir, HEADCROP=0, logger=False):
        """Build a command-line for trimmomatic"""
        if logger:
            if int(HEADCROP) > 0:
                self._logger.info("left crop reads at %d bases" %
                                  int(HEADCROP))
        # trimmo can trim the start of reads. Default this will be 0
        # but this can be defined by the user. 
        HEADCROP = "HEADCROP:%s" %(str(HEADCROP))
        prefix = L_reads.split("_R")[0]
        prefix = prefix.split("/")[-1]
        self._outdirname = os.path.join(outdir)
        self._Left_outfile = os.path.join(outdir,
                                          prefix + "_paired_R1.fq.gz")
        self._Right_outfile = os.path.join(outdir,
                                           prefix + "_paired_R2.fq.gz")

        cmd = ["java", "-jar",
               trimmo_prog,
               "PE",
               "-threads", str(threads),
               "-phred33",
               L_reads, R_reads,
               self._Left_outfile,
               "unpaired_R1.fq.gz",
               self._Right_outfile,
               "unpaired_R2.fq.gz",
               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3",
               HEADCROP,
               "TRAILING:3", "SLIDINGWINDOW:4:25", "MINLEN:70"]
        self._cmd = ' '.join(cmd)


    def run(self, trimmo_prog, L_reads, R_reads, threads,
            outdir, HEADCROP="0", logger=None):
        """Run trimmomatic on the passed read files"""
        assert L_reads != R_reads, """Oh no,
        I am trying to perform trimming on two files that are the same!
        Something has gone wrong in determining the left and right
        read files."""
        self.__build_cmd(trimmo_prog, L_reads, R_reads, threads,
                         outdir, HEADCROP)
        if not os.path.exists(self._outdirname):
            if logger:
                self._logger.info("Creating output directory: %s" %
                                  self._outdirname)
            os.makedirs(self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        # for m in msg:
        #    self._logger.info(m)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if pipe.returncode != 0:
            if logger:
                self._logger.error("trimmomatic generated some errors")
            sys.exit(1)
        if pipe.returncode == 0:
            if logger:
                self._logger.info(pipe)
        print ("pipe.args = ", pipe.args)
        return (self._outdirname, pipe.args)
