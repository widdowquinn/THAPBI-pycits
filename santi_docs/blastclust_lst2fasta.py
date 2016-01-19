#!/usr/bin/env python
from Bio import SeqIO


contador=0

for linea in open("fastqjoin.join.blastclust99.lst"):
    lista_linea=[]
    contador=contador+1
    print contador
    lista=[]
    linea=linea.split()
    if len(linea)>1:    
        lista.extend(linea)
        for seq_read in SeqIO.parse(open("fastqjoin.join.fasta"),"fasta"):
            if seq_read.id in lista:
                print seq_read.id

                lista_linea.append(seq_read)

    SeqIO.write(lista_linea,open("fastqjoin.join.blastclust99_OTU"+str(contador)+".fasta","w"),"fasta")
    
