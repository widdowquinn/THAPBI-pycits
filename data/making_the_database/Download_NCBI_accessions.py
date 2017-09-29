from Bio import Entrez, SeqIO
import pandas as pd
import re

# Set your email
Entrez.email = "peter.thorpe@hutton.ac.uk"

# Load the ITS excel spreadsheet
its_data = pd.read_excel('species_list_080615_updatedMar2017_forPete.xlsx')

# Identify the GenBank accessions and restrict to
# correctly-formatted accessions
accessions = [print(acc.strip()) for acc in its_data['ITS GenBank Accession No'] if re.match('[A-Z]{2}[0-9]{6}', str(acc.rstrip()))]
print("Properly-formatted accessions: %d" % len(accessions))

# Acquire the NCBI records for each accession as SeqRecord objects
idlist = ",".join(accessions)
handle = Entrez.efetch(db="nucleotide", id=idlist, rettype="gb",
                       retmode="text")
records = list(SeqIO.parse(handle, "genbank"))
print("Downloaded %d GenBank records from NCBI" % len(records))

for record in records:
    print("Accession: {0}, sequence length: {1}".format(record.id,
                                                        len(record.seq)))
