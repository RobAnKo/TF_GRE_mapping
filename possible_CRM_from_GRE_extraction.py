#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 16:26:13 2018

@author: rkoch
"""

import json

GRE_CRMs = {}
for i in range(1, 3546):
    GRE_CRMs[i] = json.load(open('out%s.txt'%i, 'r'))


#loop over all GREs
for i in GRE_CRMs.keys():
    CRMs = GRE_CRMs[i]
    no_of_CRMs = len(CRMs)
    #loop over all CRMs in the GRE
    for j in range(no_of_CRMs):
        CRM = CRMs[j]
        PSSM = CRM['pssm']
        motif_len = len(PSSM)