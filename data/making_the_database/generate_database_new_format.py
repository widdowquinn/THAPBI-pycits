#!/usr/bin/env python
# script to return the sequences from a tab separated table of
# iTS enteries of interest.
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe
import os
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
from optparse import OptionParser



def get_table_data(infile):
    """function to open a tab sep table and return
    info required.
    Returns dictionaries.
    species_to_clade
    accession_to_species"""
    accession_to_species = dict()
    species_to_clade = dict()
    accession_list = []
    with open(infile) as handle:
        for line in handle:
            # print(line)
            if line.startswith("#"):
                continue
            if not line.strip():
                continue # if the last line is blank
            Date_described, Species_name, ITS_Clade, \
                            sub_clade, Country, \
                            Geographical_distribution, \
                            Recorded_Host_Substrate, \
                            Preferred_environment, \
                            DNA_only, Comment, Temp_Optima, \
                            Temp_Maxima, Comments, \
                            ITS_GenBank_Accession_No,\
                            Cox, Additonal\
                            = line.split("\t")
            Species_name = Species_name.replace(" ", "_")
            accession_to_species[ITS_GenBank_Accession_No.rstrip()] = Species_name.rstrip()
            clade_sub_clade = ITS_Clade + sub_clade
            species_to_clade[Species_name.rstrip()] = clade_sub_clade.rstrip()
            accession_list.append(ITS_GenBank_Accession_No.rstrip())
    return accession_to_species, species_to_clade, accession_list


def harvest_desciption_for_info(descip):
    """function to split to line up at extra > sympbols.
    Returns a list split at the > symbol"""
    return descip.split(">")


def get_info(fasta_accession, accession_to_species, species_to_clade):
    """func to return info required. Takes in tow dictionaries"""
    species = accession_to_species[fasta_accession.split(".")[0]]
    clade = species_to_clade[species]
    seq_id = "%s_%s_%s" % (clade,
                           species,
                           fasta_accession.split(".")[0])
    discription = ""
    return species, clade, seq_id, discription


def seq_getter(dbfasta, tab_file, out_file):
    "fnction to count total exome size"
    accession_to_species, species_to_clade, accession_list\
                          = get_table_data(tab_file)
    f_out = open(out_file, 'w')
    wanted_set = set([])
    for entry in accession_list:
        wanted_set.add(entry.split(".")[0])
    name_set = set([])
    need = open("all_accessions.txt", "w")
    for seq_record in SeqIO.parse(dbfasta, "fasta"):
        fasta_accession = seq_record.id.split("|")[3]
        need.write(fasta_accession + "\n")
        if fasta_accession.split(".")[0] in wanted_set:
            name_set.add(fasta_accession.split(".")[0])
            species, clade, seq_record.id, \
                     seq_record.description = get_info(fasta_accession,
                                                       accession_to_species,
                                                       species_to_clade)
            SeqIO.write(seq_record, f_out, "fasta")
        # NOT else, and allow to continue. Names are inbeded in the
        # Description of the record.
        elements = harvest_desciption_for_info(seq_record.description)
        for ele in elements:
            if ele == "": continue
            fasta_accession = ele.split("|")[3]
            if fasta_accession.split(".")[0] in wanted_set:
                name_set.add(fasta_accession.split(".")[0])
                species, clade, seq_record.id, \
                     seq_record.description = get_info(fasta_accession,
                                                       accession_to_species,
                                                       species_to_clade)
                SeqIO.write(seq_record, f_out, "fasta")

    not_found = wanted_set.difference(name_set)
    print("Wanted number = %d" % (len(wanted_set)))
    print("we found %d" % (len(wanted_set.intersection(name_set))))
    print("did not get %s" % not_found)
    print("\n")
    print("WANTED:")
    for missing in not_found:
        species = accession_to_species[missing]
        clade = species_to_clade[species]
        print("%s_%s_%s" % (clade, species, missing))
    f_out.close()


usage = """Use as follows:
$ generate_database.py -i in.fasta -t tabfile.txt -o out.fasta
requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-t",  dest="tab_file", default=None,
                  help="tab_file containing the database informtation")
parser.add_option("-d",  dest="dbfasta", default=None,
                  help="database fasta file to get seq from")


(options, args) = parser.parse_args()

out_file = options.out_file
tab_file = options.tab_file
dbfasta = options.dbfasta

(options, args) = parser.parse_args()

if __name__ == '__main__':
    seq_getter(dbfasta, tab_file, out_file)


