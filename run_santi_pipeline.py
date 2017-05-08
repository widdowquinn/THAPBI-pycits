#!/usr/bin/env python3
#
# run_santi_pipeline.py
#
# Script to identify OTUs from metabarcoding reads.
#
# This is an almost direct translation of a pipeline written by Santiago
# Garcia (Forest Research), to generate OTU clusters from metabarcoding ITS
# reads in Phytophthora.
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe

import sys

from pycits.scripts import santi_pipeline

# Run as script
if __name__ == '__main__':

    sys.exit(santi_pipeline.run_pycits_main())
