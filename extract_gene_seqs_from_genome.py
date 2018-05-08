#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 17:41:52 2018

@author: rkoch
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:39:18 2018

@author: rkoch


Here we extract all sequences for all genes



"""
import os
from Bio import SeqIO
import pickle
os.chdir('/Users/rkoch/Documents/Data_and_Scripts')

with open("MTB_genome_data_12/genome.fna") as f:
    record = SeqIO.parse(f, "fasta")
    genome = list(record)[0]


with open('./Data/gene_data.csv', 'r') as f:
    lines = f.readlines() 
gene_positions = {}                                                          #gene: boundaries
for line in lines[1:]:
    line = line.split(',')
    gene = line[0]
    start = int(line[1])
    end = int(line[2])
    if not gene.startswith('Rvn'):
        gene_positions[gene] = [start, end]   
    

#with positions in name
gene_seqs = [nan]*len(gene_positions)
i = 0
for gene, positions in gene_positions.items():
    start = positions[0]
    end = positions[1]+1
    sequence = genome[start:end]
    sequence.id = '%s:%s-%s'%(gene, start, end)
    sequence.name = ''
    sequence.description = ''
    gene_seqs[i] = sequence
    i+=1

gene_seqs.sort(key = lambda x: x.id)    
SeqIO.write(gene_seqs, open('MTB_genome_data_12/genes.fna', 'w'), "fasta")


###---------------------------------------Part-2-------------------------------------------------###

#a dictionary containing for every corem all genes in the corem
corem_genes = pickle.load(open('./Dictionaries/corems_genes.p', 'rb'))
#flatten und uniquify the list of genes
vs = list(corem_genes.values())
vs = list(set([gene for l in vs for gene in l]))
corem_gene_prom_seqs = []
for gene in vs:
    positions = gene_positions[gene]
    start = positions[0]
    end = positions[1]+1
    sequence = genome[start:end]
    sequence.id = '%s:%s-%s'%(gene, start, end)
    sequence.name = ''
    sequence.description = ''
    corem_gene_prom_seqs.append(sequence)
corem_gene_prom_seqs.sort(key = lambda x: x.id)
SeqIO.write(corem_gene_prom_seqs, open('MTB_genome_data_12/corem_genes.fna', 'w'), "fasta")