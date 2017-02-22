#!/usr/bin/env python
#
# Vsearch * (cluster assembled reads with database)
# follow this link to get the download.
# https://github.com/torognes/vsearch
# https://insidedna.me/tool_page_assets/pdf_manual/vsearch.pdf
# http://luckylion.de/2016/06/21/quality-filtering-metabarcoding-datasets/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
from collections import namedtuple
from .tools import is_exe, NotExecutableError

# factory class for Vsearch class returned values
Results_derep = namedtuple("Results", "command fasta " +
                           "stdout stderr")

Results_cluster = namedtuple("Results", "command blast6 uc_clusters " +
                             "stdout stderr")

Results_fasta = namedtuple("Results", "command blast6 uc_clusters " +
                           "aligned centroids consensus_cls " +
                           "stdout stderr")


class Vsearch_Error(Exception):
    """Exception raised when Vsearch fails"""
    def __init__(self, message):
        self.message = message


class Vsearch_derep(object):
    """Class for working with Vsearch dereplicate"""

    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, fasta_in, outdir, prefix, dry_run=False):
        """Run Vsearch to dereplicate the passed fasta files
        --sizeout adds the abundance and orders them in size
        order
        e.g.
        vsearch --derep_fulllength in.fasta --output out.fasta --sizeout

        Returns a tuple of output filenames, and the STOUT returned by the
        Vsearch run.
        """
        self.__build_cmd(fasta_in, outdir, prefix)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

        results = Results_derep(self._cmd, self._outfnames,
                                pipe.stdout,
                                pipe.stderr)
        return results

    def __build_cmd(self, fasta_in, outdir, prefix):
        """Build a command-line for Vsearch dereplicate.

        Vsearch takes an input fasta
        and an ouput  - outdir

        path to an output directory PLUS the prefix of the
        files to write, such that

        -o a/b/cdefg
        """
        # outfiles are name WhatEver.out + .bak.clstr and + .clstr
        self._outfnames = os.path.join(outdir, prefix +
                                       'derep.fasta')
        cmd = ["vsearch",
               "--derep_fulllength", fasta_in,
               "--output", self._outfnames,
               "--sizeout"]
        self._cmd = ' '.join(cmd)

#######################################################################


class Vsearch_cluster(object):
    """Class for working with Vsearch cluster"""

    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, fasta_in, outdir, prefix, db,
            threads, threshold=0.99, dry_run=False):
        """Run Vsearch to cluster the passed fasta files
        - the fasta file should already have been dereplicated and
        sorted using the above class.


        Returns a tuple of output filenames, and the STOUT returned by the
        Vsearch run.
        """
        self.__build_cmd(fasta_in, outdir, prefix, db, threads, threshold)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

        results_clu = Results_cluster(self._cmd, *self._outfnames,
                                      pipe.stdout,
                                      pipe.stderr)
        return results_clu

    def __build_cmd(self,  fasta_in, outdir, prefix, db, threads,
                    threshold):
        """Build a command-line for Vsearch_cluster.

        Vsearch takes an input fasta
        and an ouput  - outdir

        path to an output directory PLUS the prefix of the
        files to write, such that

        -o a/b/cdefg
        """
        threshold = str(threshold)
        self._outfnames = [os.path.join(outdir, prefix) + suffix for suffix in
                           (threshold + '.blast6',
                            threshold + '.clusters.uc')]

        cmd = ["vsearch",
               "--usearch_global",
               fasta_in,
               "--id",
               threshold,
               "--uc",
               os.path.join(outdir, prefix + threshold + '.clusters.uc'),
               "--db",
               db,
               "--threads",
               str(threads),
               "--blast6out",
               os.path.join(outdir, prefix + threshold + '.blast6')]
        self._cmd = ' '.join(cmd)


class Vsearch_fastas(object):
    """Class for working with Vsearch extra tools"""

    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, fasta_in, outdir, prefix, db,
            threads, threshold=0.99, dry_run=False):
        """Run Vsearch to cluster the passed fasta files
        - the fasta file should already have been dereplicated and
        sorted using the above class.
        this specific extra class runs an alignment and return
        representative seq for cluster too


        Returns a tuple of output filenames, and the STOUT returned by the
        Vsearch run.
        """
        self.__build_cmd(fasta_in, outdir, prefix, db, threads, threshold)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

        results_fas = Results_fasta(self._cmd, *self._outfnames,
                                    pipe.stdout,
                                    pipe.stderr)
        return results_fas

    def __build_cmd(self,  fasta_in, outdir, prefix, db, threads,
                    threshold):
        """Build a command-line for Vsearch_fastas.

        Vsearch takes an input fasta
        and an ouput  - outdir

        path to an output directory PLUS the prefix of the
        files to write, such that

        -o a/b/cdefg
        """
        self._outfnames = [os.path.join(outdir, prefix + threshold) +
                           suffix for suffix in
                           ('.clusterfast.blast6',
                            '.fast.clusters.uc',
                            '.alignedclusters.fasta',
                            '.centroids.fasta',
                            '.consensus_cls_seq.fasta')]

        cmd = ["vsearch",
               "--cluster_fast",
               fasta_in,
               "--id",
               str(threshold),
               "--centroids",
               os.path.join(outdir, prefix + threshold +
                            '.centroids.fasta'),
               "--msaout",
               os.path.join(outdir, prefix + threshold +
                            '.alignedclusters.fasta'),
               "--uc",
               os.path.join(outdir, prefix + threshold +
                            '.fast.clusters.uc'),
               "--consout",
               os.path.join(outdir, prefix + threshold +
                            '.consensus_cls_seq.fasta'),
               "--db",
               db,
               "--threads",
               str(threads),
               "--blast6out",
               os.path.join(outdir, prefix + threshold +
                            '.clusterfast.blast6')]
        self._cmd = ' '.join(cmd)
