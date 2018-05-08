#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 19:18:34 2018

@author: rkoch


analysis of motif-based matching
"""

import os
import pickle
import pandas as pd
nan = float('NaN')
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages

combined_match_list_all = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/combined_match_list_all.p', 'br'))
combined_match_list_in = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/combined_match_list_in.p', 'br'))
combined_match_list_out = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/combined_match_list_out.p', 'br'))
combined_match_list_in_and_expressed = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/combined_match_list_in_and_expressed.p', 'br'))

lists = {'in':combined_match_list_in, 'in_ex':combined_match_list_in_and_expressed, 'out':combined_match_list_out, 'all': combined_match_list_all}

"""
#make unique IDs based on the TF ID, GRE ID and motif ID 
for k in lists.keys():
    lists[k] = [['^'.join([entry[0], entry[1], str(entry[2])]), entry[3], entry[4], entry[5], entry[6], entry[7], entry[8]] for entry in lists[k]]
"""    

#store all the results as dataframe
for k in lists.keys():
    lists[k] = pd.DataFrame(lists[k], columns = ['GRE_ID','TF_ID', 'TF shape', 'TF motif ID','case_overlap', 'length_overlap', 'nts_cut', 'pn', 'pc', 'orientation'])

#create match lists based on pvalue threshold 
sign_lists = {}
for k in lists.keys():
    sig = lists[k].pc < 0.05
    sign_lists[k] = lists[k][sig]

#add ranks for all values by forst sorting, then adding the ranks with the indices of sorted df
for key in lists.keys():
    lists[key]['ratio'] =  lists[key].pn/lists[key].pc     
    lists[key] = lists[key].sort_values(by = "pn")
    lists[key]['pn_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pn_rank', index = lists[key].index)  
    lists[key] = lists[key].sort_values(by = 'pc')
    lists[key]['pc_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pc_rank', index = lists[key].index)  
    lists[key] = lists[key].sort_values(by = 'ratio', ascending = False)
    lists[key]['ratio_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'ratio_rank', index = lists[key].index)
    #also for filtered list
    sign_lists[key]['ratio'] =  sign_lists[key].pn/sign_lists[key].pc     
    sign_lists[key] = sign_lists[key].sort_values(by = "pn")
    sign_lists[key]['pn_rank'] = pd.Series(range(1, len(sign_lists[key])+1), name = 'pn_rank', index = sign_lists[key].index)  
    sign_lists[key] = sign_lists[key].sort_values(by = 'pc')
    sign_lists[key]['pc_rank'] = pd.Series(range(1, len(sign_lists[key])+1), name = 'pc_rank', index = sign_lists[key].index)  
    sign_lists[key] = sign_lists[key].sort_values(by = 'ratio', ascending = False)
    sign_lists[key]['ratio_rank'] = pd.Series(range(1, len(sign_lists[key])+1), name = 'ratio_rank', index = sign_lists[key].index)
    
    

#plot pn ranks vs pc ranks
pp = PdfPages('filtered_ranks_pn_vs_pc.pdf')   
for k in sign_lists.keys():
    f = pl.figure()
    ax = sign_lists[k].plot.scatter('pn_rank', 'pc_rank',s= 0.001, title = 'ranks of normal vs cut pvalues in dataset %s'%k)
    pp.savefig(ax.figure)

pp.close()

pp = PdfPages('filtered_pn_vs_pc.pdf')   
for k in sign_lists.keys():
    f = pl.figure()
    ax = sign_lists[k].plot.scatter('pn', 'pc',s= 0.001, title = 'normal vs cut pvalues in dataset %s'%k, logx = True, logy = True)
    pp.savefig(ax.figure)

pp.close()

pp = PdfPages('filtered_ratio_rank_vs_pn.pdf')   
for k in sign_lists.keys():
    f = pl.figure()
    ax = sign_lists[k].plot.scatter('ratio_rank', 'pn',s= 0.001, title = 'ratio rank vs normal pvalues in dataset %s'%k, logx = False, logy = True)
    pp.savefig(ax.figure)

pp.close()


    


##ranking for all sorted by normal pvalue or cut pvalue
#ranks_a_n = [nan]* len(all_sorted_n)
#for i in range(len(all_sorted_n)):
#    entry = all_sorted_n[i]
#    ranks_a_n[i] = [''.join([entry[0], entry[1], str(entry[2])]), i]
#    
#ranks_a_c = [nan]* len(all_sorted_c)
#for i in range(len(all_sorted_c)):
#    entry = all_sorted_c[i]
#    ranks_a_c[i] = [''.join([entry[0], entry[1], str(entry[2])]), i]
#
#    
#
##ranking for in and expressed sorted by normal pvalue or cut pvalue
#ranks_ie_n = [nan]* len(in_ex_sorted_n)
#for i in range(len(in_ex_sorted_n)):
#    entry = in_ex_sorted_n[i]
#    ranks_ie_n[i] = [''.join([entry[0], entry[1], str(entry[2])]), i]
#    
#ranks_ie_c = [nan]* len(in_ex_sorted_c)
#for i in range(len(in_ex_sorted_c)):
#    entry = in_ex_sorted_c[i]
#    ranks_ie_c[i] = [''.join([entry[0], entry[1], str(entry[2])]), i]





##store naive ranks for all four sets for every match
#ranks = []
#indexi = 0
#for i in [combined_match_list_all, combined_match_list_in, combined_match_list_in_and_expressed, combined_match_list_out]:
#    ranks[indexi] = [nan]*len(i)
#    indexj = 0
#    for j in i:
#        ID = ''.join(j[0:3])
#        ranks[indexi][indexj] = [ID, indexj]
#        indexj += 1
#    indexi += 1
 