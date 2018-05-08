#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 18:29:03 2018

@author: rkoch
"""

reciproc_matches = {}
for key in sublists_TF.keys():
    reciproc_matches[key] = []
    progling = pb.ProgressBar()
    for n in progling(range(1,10)):
        top_TFs = {}
        top_GREs = {}
        for GRE in GRE_IDs:
            top_TFs[GRE] = sublists_GRE[key][GRE].iloc[0:n]
        for TF in TF_IDs:
            top_GREs[TF] = sublists_TF[key][TF].iloc[0:n]
        for G, Ts in top_TFs.items():
            for T in Ts.TF_ID:
                if G in top_GREs[T].GRE_ID.values:
                    reciproc_matches[key].append((G,T,n))
                    print('yaay')



no_reciprocal_matches = {}
for key in sublists_TF.keys():
    no_reciprocal_matches[key] = [len(x) for x in reciproc_matches[key]]


unique_reciprocal_matches = set([x for k in reciproc_matches.keys() for n in reciproc_matches[k].keys() for x in reciproc_matches[k][n]])
reciproc_GREs = [x[0] for x in unique_reciprocal_matches]
reciproc_TFs = [x[1] for x in unique_reciprocal_matches]