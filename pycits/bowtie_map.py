#!/usr/bin/env python
#
# Tools for working with Bowtie2 map
#
# http://bowtie-bio.sourceforge.net/bowtie2/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
from collections import namedtuple
from .tools import is_exe, NotExecutableError

# factory class for bowtie build class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# command -  the command used for the mapping
# sam - this is the outputed sam file.
# stderr
Results = namedtuple("Results", "command sam stdout stderr")


class Bowtie2_Map(Exception):
    """Exception raised when Bowtie2_Map fails"""
    def __init__(self, message):
        self.message = message


class Bowtie2_Map(object):
    """Class for working with Bowtie2_Map"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, reads, index_prefix, outdir, threads,
            dry_run=False):
        """
        Bowtie2 is btter for longer reads, than bowtie 1.
        Run Bowtie2_Map on the passed file
        reads       - can be a fasta file or a string of left and right reads.
        infnames    -  the fasta to index
        index_out   -  the out index prefix

        if perfect mapping is required consider:
        --score-min 'C,0,-1'
        """
        self.__build_cmd(reads, index_prefix, outdir, threads)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfname, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, reads, index_prefix, outdir, threads):
        """Build a command-line for Bowtie2_Map
        Bowtie2 take many option.
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#
        For the mapping we will ask for very sensitive and
        filter the sam output for the number of mismatches
        we are interested in.

        --very-sensitive    -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
        Some forums say that -N is number of mismatches. THIS IS WRONG

        --no-unal    do not report unaligned reads. Keep the size of the file
        down.

        -p     Threads
        -x     this is the index prefix which must be set up before this is
        run

        -f/-r    -f is a fasta file, -r are comma separated read lists.
        The script currently guesses if it is a fasta ot fq file.

        -S     this is the output sam file. script return the sam file as
        the fast or fq name upto the file extention and the prefix for the
        databse searched.
        eg.

        -x targets -f million_pound_gene_db.fasta -S would be ...
         -S million_pound_gene_db_Vs_targets.sam

         We will filter the output on the number of mismatches:
         XM:i:<n> The number of mismatches in the alignment

        """
        first = os.path.split(reads)[-1].split(".f")[0]
        second = "_Vs_" + os.path.split(index_prefix)[-1] + ".sam"
        outfile_name = first + second
        self._outfname = os.path.join(outdir, outfile_name)
        file_type = ""
        for i in reads.split():
            # if the reads are fasta or fasta.gz need -f
            if i.replace(".gz", "").endswith(".fa") or \
               i.replace(".gz", "").endswith(".fasta"):
                # if fasta file
                file_type = "-f"
            else:
                # else these are reads -r
                file_type = "-q"
        cmd = ["bowtie2",
               "--very-sensitive",
               "--no-unal",
               "-p", threads,
               "-x",
               index_prefix,
               file_type,
               reads,
               "-S",
               self._outfname]
        self._cmd = ' '.join(cmd)
