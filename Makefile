# Makefile
#
# Reimplementation of Santiago Garcia's OTU clustering pipeline in Make
#
# Input at the command line would be by
#
# make -e PREFIX=<path to prefix of sequence data> SEQDIR=<path to dir>
#
# We're assuming paired-end reads, where the PREFIX ends with
# - _R1_001.fastq.gz
# - _R2_001.fastq.gz
# for forward and reverse reads, respectively
#
# This analysis assumes we want to run one prefix sequence data set at a time.
# There are several third-party dependencies, and dependence on two 
# Python files that come from Santi.
# 
# Author: Leighton Pritchard
# (c) 2016 The James Hutton Institute

#PREFIX=DNAMIX_S95_L001
# Largest dataset:
PREFIX=SRB-150114_S61_L001
# Smallest dataset:
#PREFIX=SRB-030713_S56_L001

SEQDIR=../2015-12-14_THAPBI_data/MiSeq_Run2_gz_fastq_Results
SEQPREF=$(SEQDIR)/$(PREFIX)
INFILES=$(wildcard $(SEQPREF)*.gz)
FASTQCRAW_DIR=./fastqc_$(PREFIX)
FASTQC_RAW=$(patsubst $(SEQPREF)%.fastq.gz, $(FASTQCRAW_DIR)/$(PREFIX)%_fastqc.html, $(INFILES))
TRIMMED=$(patsubst $(SEQPREF)%.fastq.gz, $(PREFIX)%_trimmed.fastq, $(INFILES))
FASTQCTRIM_DIR=./fastqc_$(PREFIX)_trimmed
FASTQC_TRIM=$(patsubst $(PREFIX)%.fastq, $(FASTQCTRIM_DIR)/$(PREFIX)%_fastqc.html, $(TRIMMED))
FASTQCJOIN_DIR=./fastqc_$(PREFIX)_joined
FASTQJOINED=./fastqjoin.join
ALIGN_DIR=./align

# Run clustering
.PHONY : otus
otus : aligned_OTUs threshold_qiime_OTUs fastqc_raw fastqc_trimmed fastqc_joined closed_qiime_OTUs 


# Run FASTQC on raw reads for $(PREFIX)
.PHONY : fastqc_raw
fastqc_raw : $(FASTQC_RAW)

$(FASTQCRAW_DIR)/%_fastqc.html : $(SEQDIR)/%.fastq.gz
	mkdir -p $(FASTQCRAW_DIR)
	fastqc -o $(FASTQCRAW_DIR) $^


# Trim raw reads for quality and write trimmed reads to local directory
.PHONY : trim_raw
trim_raw : $(TRIMMED)

%_trimmed.fastq : $(SEQDIR)/%.fastq.gz
	trim_quality -t fastq --paired_reads -o $@ $^


# Run FASTQC on trimmed reads for $(PREFIX)
.PHONY : fastqc_trimmed
fastqc_trimmed : $(FASTQC_TRIM)

$(FASTQCTRIM_DIR)/%_fastqc.html : %.fastq
	mkdir -p $(FASTQCTRIM_DIR)
	fastqc -o $(FASTQCTRIM_DIR) $^


# Join trimmed paired-end files
$(FASTQJOINED).fastq : $(TRIMMED)
	join_paired_ends.py -f $(word 1, $^) -r $(word 2, $^) -o ./


# Run FASTQC on joined, trimmed reads for $(PREFIX)
.PHONY : fastqc_joined
fastqc_joined : $(FASTQJOINED).fastq
	mkdir -p $(FASTQCJOIN_DIR)
	fastqc -o $(FASTQCJOIN_DIR) $<


# Convert joined FASTQ to FASTA and trim with 'trim_longitudes.py'
$(FASTQJOINED).fasta : $(FASTQJOINED).fastq
	convert_format -t fastq -o $(FASTQJOINED).fasta \
		       -f fasta $(FASTQJOINED).fastq
	trim_longitudes.py


# Cluster FASTA reads with BLASTCLUST
.PHONY : cluster
cluster : $(FASTQJOINED).blastclust99.lst

%.blastclust99.lst : %.fasta
	blastclust -L 0.90 -S 99 -a 4 -p F \
	           -i $(FASTQJOINED).fasta -o $@ 


# Divide BLASTCLUST output into separate files
.PHONY : blast_OTUs
blast_OTUs : $(ALIGN_DIR)/$(FASTQJOINED).blastclust99_OTU1.fasta

$(ALIGN_DIR)/$(FASTQJOINED).blastclust99_OTU1.fasta : $(FASTQJOINED).blastclust99.lst
	mkdir -p align
	blastclust_lst2fasta.py
	mv *OTU* $(ALIGN_DIR)


# Align clusters with MUSCLE
.PHONY : aligned_OTUs
aligned_OTUs : $(ALIGN_DIR)/$(FASTQJOINED).blastclust99_OTU1.fasta.aln

$(ALIGN_DIR)/$(FASTQJOINED)%.aln : $(ALIGN_DIR)/$(FASTQJOINED)%
	for filename in $(ALIGN_DIR)/*.fasta; \
	do \
	  (muscle -in $$filename -out $$filename.aln) \
	done


# Pick OTUs from reference FASTA database with QIIME
.phony : closed_qiime_OTUs
closed_qiime_OTUs : otus/otu_table.tab

otus/otu_table.tab : fastqjoin.join.fasta reference.fasta
	pick_closed_reference_otus.py -f -i fastqjoin.join.fasta \
	                              -r reference.fasta -o otus/
	biom convert -i otus/otu_table.biom -o otus/otu_table.tab --to-tsv


# Pick OTUs with QIIME on basis of similarity threshold
.phony : threshold_qiime_OTUs
threshold_qiime_OTUs : uclust_ref_picked_otus

uclust_ref_picked_otus : fastqjoin.join.fasta reference.fasta
	pick_otus.py -i fastqjoin.join.fasta -r reference.fasta \
	             -m uclust_ref -s 0.99 --threads=8


# Clean up output
.PHONY : clean
clean:
	rm -rf fastqc*
	rm *_trimmed.fastq
	rm fastqjoin*
	rm -rf align
	rm -rf otus
	rm -rf uclust_ref_picked_otus

