#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:55:29 2018
Due to the idea to run tomtom not against every single GRE but against a database of all GREs
@author: rkoch
"""

import os
import re



os.chdir('/Users/rkoch/Documents/Data_and_Scripts/motif_clusters_24/')

files = os.listdir()
files = [x for x in files if x.endswith('memeOut.txt')]
merged_file = open('merged_GRE_motifs.txt', 'w')
first_flag = True
for file in files:
    ID = int(file.split('_')[0]) + 1
    with open(file) as f:
        lines = f.readlines()
        line1 = [x for x in lines if x.startswith('MOTIF  1 MEME	width')][0]    # header line
        i1 = lines.index(line1)                                                    # index of header line
        n = int(line1.split()[5])                                                  # size of motif
        line1 = re.sub('\d+', str(ID), line1,  count =1)                                                  
        block1 = [lines[i1-1], line1, lines[i1+1]]                                 # first block to introduce into big file
        
        line2 = [x for x in lines if x.startswith('log-odds matrix:')][0]
        i2 = lines.index(line2)
        block2 = lines[i2-1:i2+n+1]
        
        if first_flag:
            head = lines[0:4]
            for i in head:
                merged_file.write(i)
            first_flag = False
            
        for i in block1:
            merged_file.write(i)
        for i in block2:
            merged_file.write(i)
merged_file.close()