READme Making and clustering an ITS1 database
=========================
The database from Santii's oringinal pipeline was aligned with a known ITS1 region from P. infestans obtained from NCBI GenBank.
The coordintaes of the ITS1 regions were defined in GenBank. The original database had sequences of around 800nt, 
which would not be suitable for clustering with swarm. 
Therefore, following alignmnet with Muscle (with the --refine option), the alignment was trimmed to only the ITS1 region for all sequences. 
Thus yielding a database of ITS1 regions only.

The sequences start (with some variation) "CCACACC" and end with some variation "TTGC"