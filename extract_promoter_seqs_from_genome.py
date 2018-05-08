#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:39:18 2018

@author: rkoch


Here we extract all 220 bp promoter sequences for all genes. Part 2 does the same only for genes in corems



"""
import os
from Bio import SeqIO
import pickle
os.chdir('/Users/rkoch/Documents/Data_and_Scripts')

with open("MTB_genome_data_12/83332.12.fna") as f:
    record = SeqIO.parse(f, "fasta")
    genome = list(record)[0]


gene_dic = pickle.load(open('Dictionaries/gene_dic.p', 'br'))

"""
#without positions in name
gene_prom_seqs = []

for gene, positions in gene_dic.items():
    sequence = genome[positions[0]:positions[1]+1]
    sequence.id = gene
    sequence.name = ''
    sequence.description = ''
    gene_prom_seqs.append(sequence)

gene_prom_seqs.sort(key = lambda x: x.id)    
SeqIO.write(gene_prom_seqs, open('MTB_genome_data_12/gene_promotors.fna', 'w'), "fasta")
"""       

#with positions in name
gene_prom_seqs = []

for gene, positions in gene_dic.items():
    start = positions[0]
    end = positions[1]+1
    sequence = genome[start:end]
    sequence.id = '%s:%s-%s'%(gene, start, end)
    sequence.name = ''
    sequence.description = ''
    gene_prom_seqs.append(sequence)

gene_prom_seqs.sort(key = lambda x: x.id)    
SeqIO.write(gene_prom_seqs, open('MTB_genome_data_12/gene_promotors.fna', 'w'), "fasta")


###---------------------------------------Part-2-------------------------------------------------###

#a dictionary containing for every corem all genes in the corem
corem_genes = pickle.load(open('./Dictionaries/corems_genes.p', 'rb'))
#flatten und uniquify the list of genes
vs = list(corem_genes.values())
vs = list(set([gene for l in vs for gene in l]))
corem_gene_prom_seqs = []
for gene in vs:
    positions = gene_dic[gene]
    start = positions[0]
    end = positions[1]+1
    sequence = genome[start:end]
    sequence.id = '%s:%s-%s'%(gene, start, end)
    sequence.name = ''
    sequence.description = ''
    corem_gene_prom_seqs.append(sequence)
corem_gene_prom_seqs.sort(key = lambda x: x.id)
SeqIO.write(corem_gene_prom_seqs, open('MTB_genome_data_12/corem_gene_promotors.fna', 'w'), "fasta")