READme Making and clustering an ITS1 database
=========================
The database from Santii's oringinal pipeline was aligned with a known ITS1 region from P. infestans obtained from NCBI GenBank.
The coordintaes of the ITS1 regions were defined in GenBank. The original database had sequences of around 800nt, 
which would not be suitable for clustering with swarm. 
Therefore, following alignmnet with Muscle (with the --refine option), the alignment was trimmed to only the ITS1 region for all sequences. 
Thus yielding a database of ITS1 regions only.

The sequences start (with some variation) "CCACACC" and end with some variation "TTGC"


Problems:
>7Phytophthora_melonis_CBS58269
CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACAAGTAGTTGGGGGCCTGCTCT
GTGTGGCTAGCTGTCGATGTCAAAGTCGGCGACTGGCTGCTATGTGGCGG
GCTCTATCATGGCGATTGGTTTGGGTCCTCCTCGTGGGGAACTGGATCATGAG
CCCACCTTTTAAACCCATTCTTGATTACTGAATATACTGTGGGGACGAAAGTCTCTGC
>7Phytophthora_sinensis_CBS55788

These have the same sequence and are now named as: >7Phytophthora_melonis_CBS58269_7Phytophthora_sinensis_CBS55788

Problem 2: Duplicates have been removed.

Problem3: 
8Phytophthora_austrocedri is now represented by this sequence rather than 5, some of which were duplicates.
>8Phytophthora_austrocedri_TDJ3_ITS                   
CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAATCCTTTTATTGGGGGCTTCTGT
CTGGTCTGGCTTCGGCTGGTCTGGGTGGCGGCTCTATCATGGTGACCGCTCTGGGCTTCG
GCTTGGAGTTAGTAGCCCACTTTTTAAACCCATTCTTAATTACTGAACATACTGTGGGGA
CGAAAGTCTCTGC

Problem 4:
duplicate seq in different species again. and renamed as:
>3Phytophthora_quercina_CBS78195_3Phytophthora_sp_ohioensis_P16050
CCACACCTAAAAAAACTTTCCACGTGAACCGTTTCAACCAAATATTTTGGGGGTCTTGTC
TGGCGTATGGCTGCTGCTGTAAAAGGCGGCGGCTGTTGCTGGGTGA
GCCCTATCATGGCGAACGTTTGGGCTTCGGTCTGAACAAGTAG
CTCTTTTTTAAACCATTACTTATTACTGATTATACTGTGGGGACGAAAGTCTCTGC
>3Phytophthora_sp_ohioensis_P16050

problem 5:
Duplicate sequence found with ID
8Phytophthora_himalayensis_CBS35759

>8Phytophthora_erythroseptica_CBS111343
CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCTTTTTAAATTGGGGGCTTCCGTC
TGGCCGGCCGGTTTTCGGCTGGCTGGGTGGCG
GCTCTATCATGGCGACCGCTTGGGCCTCGGCCTGGGCTAGTAG
CGTATTTTTAAACCATTCCTAATTACTGAATATACTGTGGGGACGAAAGTCTCTGC
>8Phytophthora_himalayensis_CBS35759_8Phytophthora_erythroseptica_CBS111343
CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCTTTTTAAATTGGGGGCTTCCGTC
TGGCCGGCCGGTTTTCGGCTGGCTGGGTGGCG
GCTCTATCATGGCGACCGCTTGGGCCTCGGCCTGGGCTAGTAG
CGTATTTTTAAACCATTCCTAATTACTGAATATACTGTGGGGACGAAAGTCTCTGC