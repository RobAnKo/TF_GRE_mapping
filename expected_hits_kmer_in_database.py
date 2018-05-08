#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 09:03:27 2018

@author: rkoch
"""
from Bio import SeqIO






def exp_occurences_kmer_in_database(k, database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters.fna'):
    p = 0.25
    with open(database) as f:
                rec = SeqIO.parse(f, 'fasta')
                rec = list(rec) 
                n = sum([len(x) for x in rec])
    
    pk = p**k
    attempts = n-k+1
    return  pk * n 


ks = range(6,31)
es_promoter = []
es = []
for k in ks:
    es.append(exp_occurences_kmer_in_database(k, database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/genome.fna'))
    es_promoter.append(exp_occurences_kmer_in_database(k))
with PdfPages('expected_kmer_occurences.pdf') as pp:
    fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(ks, es, label = "whole genome")
    ax.plot(ks, es_promoter, label = "promoter regions only")
    ax.plot([6,30], [1,1], label= 'E = 1')
    ax.set_ylabel('log(Expected occurences)')
    ax.set_xlabel('length of motif k')
    ax.set_title('Expected occurences of motif with length k')
    ax.legend()
    plt.tight_layout()
    ax.set_yscale('log')
    pp.savefig()