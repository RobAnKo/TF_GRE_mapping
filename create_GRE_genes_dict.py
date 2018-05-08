#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 18:01:25 2018

@author: rkoch

Here we create a dictionary that 
"""
import os
import csv

os.chdir('/Users/rkoch/Documents/Data_and_Scripts/Data/gene_all_gres')



gre_gene = {}


files = [x for x in os.listdir() if x.endswith('.csv')]
for file in files:
    ID = file.split('.')[0][3:]
    gre_gene[ID] = []
    with open(file, 'r') as f:
        r = csv.reader(f)
        for row in r:
            if 'Rv' in row[0]:
                gre_gene[ID].extend(row)
            

pickle.dump(gre_gene, open('../../Dictionaries/gre_gene.p', 'wb'))