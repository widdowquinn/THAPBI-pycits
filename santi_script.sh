#!/usr/bin/env bash
#
# santi_script.sh
#
# Santi's pipeline - originally provided as a Word file - converted to Bash
# script.
#
# Original pipeline had no documentation or comments - all interpretation
# is new.
#
# Author: Leighton Pritchard
# (C) The James Hutton Institute 2016

# Create directory for output
mkdir fastq-join_joined

# Copy forward and reverse reads into file with 'standard' naming
cp *R1* DCM1.fastq
cp *R2* DCM2.fastq

# Conduct QC on collated "forward" and "reverse" reads
fastqc DCM1.fastq
fastqc DCM2.fastq

# Trim reads on quality
trim_quality -t fastq -o DCM1_trimmed.fastq --paired_reads DCM1.fastq
trim_quality -t fastq -o DCM2_trimmed.fastq --paired_reads DCM2.fastq

# Conduct QC on trimmed reads
fastqc DCM1_trimmed.fastq
fastqc DCM2_trimmed.fastq

# Move any Python files from the current directory to the output subdirectory
mv *.py fastq-join_joined

# Move untrimmed, collated read files to output subdirectory
mv DCM1.fastq fastq-join_joined
mv DCM2.fastq fastq-join_joined

# Join trimmed paired-end read files together
# creates file fastqjoin.join.fastq?
# Requires fastq-join (ea-utils, install via brew)
join_paired_ends.py -f DCM1_trimmed.fastq -r DCM2_trimmed.fastq \
                    -o fastq-join_joined

# Move trimmed read files and joined read files to the output subdirectory
mv *trimmed.fastq fastq-join_joined

# Remove all .zip files in current directory
rm *.zip

# Change directory to the output subdirectory
cd fastq-join_joined

# Conduct QC on joined reads
fastqc fastqjoin.join.fastq

# Convert read format from FASTQ to FASTA
# convert_format comes from seq_crumbs: https://github.com/JoseBlanca/seq_crumbs
convert_format -t fastq -o fastqjoin.join.fasta -f fasta fastqjoin.join.fastq

# Trim something?
# This seems to be adaptor sequences, from the Materials and Methods
python trim_longitudes.py

# Cluster FASTA reads with BLASTCLUST and [convert output to FASTA]?
blastclust -L 0.90 -S 99 -i fastqjoin.join.fasta \
           -o fastqjoin.join.blastclust99.lst -a 4 -p F
python blastclust_lst2fasta.py

# Create subdirectories for output
mkdir data
mkdir scripts
mv *.py scripts
mkdir align

# Move output to subdirectories
mv *OTU* align
mv *fastq data
mv *lst data
mc *.fasta data
mc *txt data
rm *.zip
mv *error* data

# Change to align subdirectory, and run MUSCLE alignment
cd align
for i in fasta
do
  muscle -in $i -out $i"_muscle.fasta"
done

# Use QIIME to pick OTUs
pick_closed_reference_otus.py -i fastqjoin.join.fasta -r reference.fasta -o otus/
pick_otus.py -i fastqjoin.join.fasta -r reference.fasta -m uclust_ref -s 0.99
