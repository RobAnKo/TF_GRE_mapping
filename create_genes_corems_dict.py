#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:48:43 2018

@author: rkoch
"""

import os
os.chdir('/Users/rkoch/egrin2-tools/src')
import egrin2.query.egrin2_query as e2q
import pymongo
import pickle

client = pymongo.MongoClient(host='como', port= 27018)

mtu_db = client['mtu_db']


corems_genes = {}
for id in range(1,561):
    genes = e2q.corem_genes(mtu_db, id)
    corems_genes[id] = genes
    
pickle.dump(corems_genes, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/corems_genes.p', 'wb'))


all_genes = []
for v in corems_genes.values():
    all_genes.extend(v)
all_genes = list(set(all_genes))

genes_corems = {}
for gene in all_genes:
    genes_corems[gene] = []
    for corem, members in corems_genes.items():
        if gene in members:
            genes_corems[gene].append(int(corem))
            
pickle.dump(genes_corems, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/genes_corems.p', 'wb'))
