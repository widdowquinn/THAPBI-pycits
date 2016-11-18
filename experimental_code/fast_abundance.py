
# coding: utf-8

# In[1]:

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import hashlib
from collections import defaultdict


# In[2]:

seqlist = list(SeqIO.parse('dedup_test.fasta', 'fasta'))
print(seqlist)


# In[3]:

for seq in seqlist:
    print(seq.id, hashlib.md5(str(seq.seq).encode()).hexdigest())


# In[4]:

hashlist = [hashlib.md5(str(seq.seq).encode()).hexdigest() for seq in seqlist]
print(hashlist)


# In[5]:

unique = set(hashlist)
print(unique, len(unique))


# In[11]:

abundance = defaultdict(int)
hash_to_seq = defaultdict(str)
hash_to_name = defaultdict(str)
for seq in seqlist:
    seqhash = hashlib.md5(str(seq.seq).encode()).hexdigest()
    abundance[seqhash] += 1
    hash_to_seq[seqhash] = seq.seq
    hash_to_name[seqhash] = seq.id
    print("{0}:\t\t\t{1}".format(seq.id, seqhash))


# In[7]:

abundance


# In[8]:

hash_to_seq


# In[12]:

for k, v in abundance.items():
    seqname = "{0}_{1}".format(k,v)
    seqrecord = SeqRecord(id=seqname, name=hash_to_name[k], seq=hash_to_seq[k])
    print(seqrecord)


# In[ ]:



