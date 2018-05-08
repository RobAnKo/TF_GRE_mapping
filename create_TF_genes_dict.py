#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:57:35 2018

@author: rkoch



Here we create a dictionary of TFs with all the genes bound by the TF in their promoter region

"""
from statistics import mean 
import pickle
import progressbar





# dictionary of all genes with values being genomic positions of their promoter region
gene_dic = pickle.load(open('Dictionaries/gene_dic.p', 'br'))


#dictionary of 156 TFs with values being a list of positions (+- 12bp around the binding position)
chip_dic = pickle.load(open('Dictionaries/chip_dic.p', 'br'))






#create a dictionary of TFs with values being the genes where the TF binds the promoter
chip_gene = {}
pbar = ProgressBar()
for TF, sites in pbar(chip_dic.items()):
    chip_gene[TF] = []
    for site in sites:
        central = mean(site)
        for gene, promoter in gene_dic.items():
            if (promoter[0] < central < promoter[1]):
                chip_gene[TF].append(gene)
                
print(chip_gene)
pickle.dump(chip_gene, open('Dictionaries/chip_genes.p', 'bw'))
                
chip_genes = pickle.load(open('Dictionaries/chip_genes.p', 'br'))
    