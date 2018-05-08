#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 14:46:15 2018

@author: rkoch
"""

import pickle
import pandas as pd
import numpy as np
import progressbar as pb
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns



#I.##--------------------------------------analysis functions--------------------------------#######


def filter_pvals(df, methods = ['uncut', 'GRE_cut', 'both_cut', 'genes', 'hits', 'nucleotides', 'peak_center'], thresh = 1):
    """
    This function takes as input a mapping DataFrame and returns a filtered DataFrame based on a pvalue threshold "thresh" for the mapping methods
    
    """
    m_dic = {'uncut': 1, 'GRE_cut': 1, 'both_cut': 1, 'genes': 1, 'hits': 1, 'nucleotides': 1, 'peak_center': 1}
    if type(methods) == list:
        for m in methods:
            m_dic[m] = thresh
    else:
        m_dic[methods] = thresh
    
    out = df.loc[(df.pv <= m_dic['uncut']) & 
                 (df.pv_cut_GRE <= m_dic['GRE_cut']) & 
                 (df.pv_cut_both <= m_dic['both_cut']) & 
                 (df.pv_genes <= m_dic['genes']) & 
                 (df.pv_hits <= m_dic['hits']) & 
                 (df.pv_nucleotides <= m_dic['nucleotides'])]
    return out


def filter_best_TF_maps(df, method = 'tomtom uncut'):
    """
    This function takes a dataframe containing all mappings and filters it for the best matches for every
    TF based on the method chosen in 'method'
    
    possible methods are:
        
        'tomtom uncut'
        'tomtom GRE cut'
        'tomtom both cut'
        'position genes'
        'position hits'
        'position nucleotides'
        'position peak center'
        'median all'
        'median motif'
        'median position'
        'best all'
        'best motif'
        'best position'
    """
    counywouny = pb.ProgressBar()
    best_TF_maps = {}
    for TF in counywouny(all_TFs):
        best_TF_maps[TF] = nan
        d = df.loc[final.TF_ID == TF]
        if method == 'tomtom uncut':
            d.sort_values(by = 'pv', inplace = True)        
        elif method == 'tomtom GRE cut':
            d.sort_values(by = 'pv_cut_GRE', inplace = True)        
        elif method == 'tomtom both cut':
            d.sort_values(by = 'pv_cut_both', inplace = True)        
        elif method == 'position genes':
            d.sort_values(by = 'pv_genes', inplace = True)        
        elif method == 'position hits':
            d.sort_values(by = 'pv_hits', inplace = True)        
        elif method == 'position nucleotides':
            d.sort_values(by = 'pv_nucleotides', inplace = True)        
        elif method == 'position peak center':
            d.sort_values(by = 'pv_peak_center', inplace = True)
        elif method == 'median all':
            d['meds'] = d[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
            d.sort_values(by = 'meds', inplace = True)   
            d.drop('meds', inplace = True, axis = 1)
        elif method == 'median motif':
            d['meds'] = d[['pv','pv_cut_GRE', 'pv_cut_both']].apply(np.median, axis = 1)
            d.sort_values(by = 'meds', inplace = True)
            d.drop('meds', inplace = True, axis = 1)
        elif method == 'median position':
            d['meds'] = d[['pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
            d.sort_values(by = 'meds', inplace = True)
            d.drop('meds', inplace = True, axis = 1)
        elif method == 'best all':
            d['best'] = d[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.min, axis = 1)
            d.sort_values(by = 'best', inplace = True)
            d.drop('best', inplace = True, axis = 1)
        elif method == 'best motif':
            d['best'] = d[['pv','pv_cut_GRE', 'pv_cut_both']].apply(np.min, axis = 1)
            d.sort_values(by = 'best', inplace = True)
            d.drop('best', inplace = True, axis = 1)
        elif method == 'best position':
            d['best'] = d[['pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.min, axis = 1)
            d.sort_values(by = 'best', inplace = True)
            d.drop('best', inplace = True, axis = 1)
        else:
            print("Not a valid method")
            return
            
        best_TF_maps[TF] = d.iloc[0]      
    out = pd.DataFrame.from_dict(best_TF_maps).transpose()
    out['index'] = out['index'].apply(int)
    out.set_index('index', inplace = True)
    return out

def filter_best_n(df, n, method = 'tomtom uncut'):
    """
    This function takes a mapping dataframe and filters for the top n hits based on the chosen method
    
    Possible methods are:
        'tomtom uncut'
        'tomtom GRE cut'
        'tomtom both cut'
        'position genes'
        'position hits'
        'position nucleotides'
        'position peak center'
        'median all'
        'median motif'
        'median position'
        'best all'
        'best motif'
        'best position'
        
    """
    #print(method)
    d = df[:]
    if method == 'tomtom uncut':
            d.sort_values(by = 'pv', inplace = True)        
    elif method == 'tomtom GRE cut':
        d.sort_values(by = 'pv_cut_GRE', inplace = True)        
    elif method == 'tomtom both cut':
        d.sort_values(by = 'pv_cut_both', inplace = True)        
    elif method == 'position genes':
        d.sort_values(by = 'pv_genes', inplace = True)        
    elif method == 'position hits':
        d.sort_values(by = 'pv_hits', inplace = True)        
    elif method == 'position nucleotides':
        d.sort_values(by = 'pv_nucleotides', inplace = True)        
    elif method == 'position peak center':
        d.sort_values(by = 'pv_peak_center', inplace = True)
    elif method == 'median all':
        d['meds'] = d[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
        d.sort_values(by = 'meds', inplace = True) 
        d.drop('meds', inplace = True, axis = 1)
    elif method == 'median motif':
        d['meds'] = d[['pv','pv_cut_GRE', 'pv_cut_both']].apply(np.median, axis = 1)
        d.sort_values(by = 'meds', inplace = True)   
        d.drop('meds', inplace = True, axis = 1)
    elif method == 'median position':
        d['meds'] = d[['pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
        d.sort_values(by = 'meds', inplace = True)
        d.drop('meds', inplace = True, axis = 1)
    elif method == 'best all':
        d['best'] = d[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.min, axis = 1)
        d.sort_values(by = 'best', inplace = True)
        d.drop('best', inplace = True, axis = 1)
    elif method == 'best motif':
        d['best'] = d[['pv','pv_cut_GRE', 'pv_cut_both']].apply(np.min, axis = 1)
        d.sort_values(by = 'best', inplace = True)    
        d.drop('best', inplace = True, axis = 1)
    elif method == 'best position':
        d['best'] = d[['pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.min, axis = 1)
        d.sort_values(by = 'best', inplace = True)    
        d.drop('best', inplace = True, axis = 1)
    else:
        print("Not a valid method")
        return
    
    return d.iloc[0:n]
    
def sig_sum(row, thresh, methods = 'both'):
    """
    a function that counts the methods for which a match (row of a dataframe) is below a certain pvalue threshold
    
    """
    
    if methods not in ['both', 'motif', 'position']:
        print("Error: methods has to be one of 'both', 'motif', 'position'")
        return
    
    ss = 0
    for n, s in row.items():
        if methods == 'both':
            if n in ['pv','pv_cut_GRE','pv_cut_both','pv_genes','pv_hits','pv_nucleotides', 'pv_peak_center']:
                if s <= thresh:
                    ss += 1
        elif methods == 'motif':
            if n in ['pv','pv_cut_GRE','pv_cut_both']:
                if s <= thresh:
                    ss += 1
        elif methods == 'position':
            if n in ['pv_genes','pv_hits','pv_nucleotides', 'pv_peak_center']:
                if s <= thresh:
                    ss += 1
    return ss
    


def GRE_genes_with_condition(condition):
    """
    get GRE-gene connections based on corems associated with condition <condition>.
    """
    GRE_genes = {}
    for GRE in unique_GREs:
        corems = GRE_corem_dic[GRE]
        for corem in corems:
            if condition != 'all':    
                if corem in corem_cond_dic:
                    if condition in corem_cond_dic[corem]:
                        if GRE in GRE_genes:
                            GRE_genes[GRE].extend(corems_genes[corem])
                        else:
                            GRE_genes[GRE] = corems_genes[corem]
            else:
                if GRE in GRE_genes:
                    GRE_genes[GRE].extend(corems_genes[corem])
                else:
                    GRE_genes[GRE] = corems_genes[corem]
    for GRE, genes in GRE_genes.items():
        GRE_genes[GRE] = set(genes)           
    return GRE_genes

#II.##-------------------------------------Import the results--------------------------------#######

data_dir = '/Users/rkoch/Documents/Data_and_Scripts/out/match_lists/'
#import the match list based on gene overlap
match_list_genes = pickle.load(open('{}match_list_genes.p'.format(data_dir), 'br'))

#import the match list based on hit overlap
match_list_hits = pickle.load(open('{}match_list_positions_hits.p'.format(data_dir), 'br'))

#import the match list based on hit nucleotide overlap
match_list_nucleotides = pickle.load(open('{}match_list_positions_nucleotides.p'.format(data_dir), 'br'))

#import the match list based on peak center overlap
match_list_peak_center = pickle.load(open('{}match_list_positions_peak_center.p'.format(data_dir), 'br'))

#transform results into DataFrames
cols = ['GRE_ID', 'TF_ID', 'n_GRE', 'n_TF', 'n_both', 'pv']
ml_g = pd.DataFrame(match_list_genes, columns = cols)
ml_g.set_index(['GRE_ID', 'TF_ID'], inplace = True)
ml_g.columns = [x + '_genes' for x in ml_g.columns]
ml_h = pd.DataFrame(match_list_hits, columns = cols)
ml_h.set_index(['GRE_ID', 'TF_ID'], inplace = True)
ml_h.columns = [x + '_hits' for x in ml_h.columns]
ml_n = pd.DataFrame(match_list_nucleotides, columns = cols)
ml_n.set_index(['GRE_ID', 'TF_ID'], inplace = True)
ml_n.columns = [x + '_nucleotides' for x in ml_n.columns]
ml_p = pd.DataFrame(match_list_peak_center, columns = cols)
ml_p.set_index(['GRE_ID', 'TF_ID'], inplace = True)
ml_p.columns = [x + '_peak_center' for x in ml_p.columns]
ml_gh = pd.merge(ml_g, ml_h, left_index = True, right_index = True)
ml_ghn = pd.merge(ml_gh, ml_n, left_index = True, right_index = True)
ml_positions = pd.merge(ml_ghn, ml_p, left_index = True, right_index = True)

#import the match lists based on motif similarity for both sequences cut
ml_all = pickle.load(open('{}combined_match_list_all.p'.format(data_dir), 'rb'))
ml_in = pickle.load(open('{}combined_match_list_in.p'.format(data_dir), 'rb'))
ml_out = pickle.load(open('{}combined_match_list_out.p'.format(data_dir), 'rb'))
ml_in_and_expressed = pickle.load(open('{}combined_match_list_in_and_expressed.p'.format(data_dir), 'rb'))

#for only GRE cut
#ml_all_s = pickle.load(open({}combined_match_list_with_center_all.p', 'rb'))
#ml_in_s = pickle.load(open({}combined_match_list_with_center_in.p', 'rb'))
#ml_out_s = pickle.load(open({}combined_match_list_with_center_out.p', 'rb'))
#ml_in_and_expressed_s = pickle.load(open({}combined_match_list_with_center_in_and_expressed.p', 'rb'))

#store all motif based match lists in a dictionary
lists = {'in':ml_in, 'in_ex':ml_in_and_expressed, 'out':ml_out, 'all': ml_all}
#short_lists = {'in':ml_in_s, 'in_ex':ml_in_and_expressed_s, 'out':ml_out_s, 'all': ml_all_s}
cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'case_of_overlap', 'length_overlap', 'length_GRE_uncut', 'length_GRE_cut', 'number_spliced_nts_GRE', 'length_TF_uncut', 'length_TF_cut', 'number_spliced_nts_TF', 'pv', 'pv_cut_both','pv_cut_GRE', 'orientation', 'center', 'evalue_TF', 'no_seqs_TF', 'evalue_GRE', 'no_seqs_GRE']
#short_cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'case_of_overlap', 'length_overlap', 'number_spliced_nts_GRE', 'length_GRE_uncut', 'length_GRE_cut', 'length_TF_uncut', 'pv', 'pv_cut', 'orientation', 'center']

#store all the results as dataframe
for k in lists.keys():
    lists[k] = pd.DataFrame(lists[k], columns = cols)
    #short_lists[k] = pd.DataFrame(short_lists[k], columns = short_cols)
    lists[k].to_csv(open('motif_matches_{}.csv'.format(k), 'w'))
    #short_lists[k].to_csv(open('motif_matches_{}_only_GRE.csv'.format(k), 'w'))


    
#III.##-----------------------------------combine the results--------------------------------#######
for key in lists.keys():
    #optional: reduce matches to those with motifID 1
    #lists[key] = lists[key].loc[lists[key].TF_motif_ID == 1] 
    #lists[key].index = range(len(lists[key]))
    
    #'improvement' of pvalue with cut GRE
    lists[key]['pv_ratio_cut_GRE'] =  lists[key].pv.divide(lists[key].pv_cut_GRE)
    
    #fraction of TF motif that is aligned with GRE motif
    lists[key]['TF_coverage'] = lists[key].length_overlap.divide(lists[key].length_TF_uncut)

    #'improvement' of pvalue with cut GRE and TF
    lists[key]['pv_ratio_cut_both'] = lists[key].pv.divide(lists[key].pv_cut_both)
    
    #add column with information w.r.t. peak set used
    lists[key]['peaks'] = len(lists[key])*[key]
    
    
#a dictionary with mappings based on motifs
motif_lists = lists.copy()

# a dataframe containing all motif-mappings
a = pd.concat(lists, ignore_index = True)
a.set_index(['GRE_ID', 'TF_ID'], inplace = True)


##-> dataframe containing all motif-mappings with position results included
complete_result = pd.merge(a, ml_positions, left_index = True, right_index = True)
#store end-result
complete_result.to_csv(open('all_matches_combined.csv','w'))
complete_result = pd.read_csv(open('all_matches_combined.csv','r'))

all_GREs = set(complete_result.GRE_ID)
all_TFs = set(complete_result.TF_ID)
all_combinations = [(x,y) for x in all_GREs for y in all_TFs]

#IV.##-----------------extract best motif-based option for each pair-------------------------#######

complete_result.reset_index(inplace = True)
df = complete_result
final1 = pd.DataFrame()
final2 = pd.DataFrame()
countyboi = pb.ProgressBar()


#one version using the minimum of uncut motif pvalues to decide on best match
for pair in countyboi(all_combinations):
    print('{} vs {}'.format(pair[0], pair[1]))
    pair_matches = df.loc[(df.GRE_ID == pair[0]) & (df.TF_ID == pair[1])]
    best_one = pair_matches.loc[pair_matches.pv == np.nanmin(pair_matches.pv)].iloc[0]
    final1 = final1.append(best_one, ignore_index = True)

countyboi = pb.ProgressBar()
#one version using the minimum of all motif pvalues to decide on best match
for pair in countyboi(all_combinations):
    print('{} vs {}'.format(pair[0], pair[1]))
    pair_matches = df.loc[(df.GRE_ID == pair[0]) & (df.TF_ID == pair[1])]
    mins = pair_matches[['pv', 'pv_cut_both', 'pv_cut_GRE']].apply(min, axis = 1)
    best_one = pair_matches.loc[mins == np.nanmin(mins)].iloc[0]
    final2 = final2.append(best_one, ignore_index = True)
 
final1.to_csv(open('best_motif_pv_mapping_combined.csv', 'w'))
final1 = pd.read_csv(open('best_motif_pv_mapping_combined.csv', 'r'))
final2.to_csv(open('best_motif_3pvs_mapping_combined.csv', 'w'))
final2 = pd.read_csv(open('best_motif_3pvs_mapping_combined.csv', 'r'))   

#writing to file and reading in does some weird stuff -> cleaning data
final1.drop('Unnamed: 0', axis = 1, inplace = True)
final1.GRE_ID = final1.GRE_ID.apply(int)
final1.TF_motif_ID = final1.TF_motif_ID.apply(int) 
final2.drop('Unnamed: 0', axis = 1, inplace = True)
final2.GRE_ID = final2.GRE_ID.apply(int)
final2.TF_motif_ID = final2.TF_motif_ID.apply(int) 

#difference between both approaches:
final11 = final1.set_index(final1.columns.tolist())
final22 = final2.set_index(final2.columns.tolist())
diff_best_pvs = final11.index.difference(final22.index)
diff_best_pvs2 = final22.index.difference(final11.index)
different_best_ones = final11.loc[diff_best_pvs] #the matches for the differing pairs based on pv
different_best_ones2 = final22.loc[diff_best_pvs2] #the matches for the differing pairs based on 3pv
different_best_ones.reset_index(inplace = True)
different_best_ones2.reset_index(inplace = True)
diff_GT = [(int(x.GRE_ID), x.TF_ID) for _,x in different_best_ones.iterrows()]

# -> preliminary result: 7800/19000 pairs have different 'best' matches between pv and 3pv approach
# are there high-confidence pairs in those 7800?
diff_GT = set(diff_GT)
t = pE1
ss_motif = final1.apply(sig_sum, axis = 1, thresh = t, methods = 'motif')
ss_position = final1.apply(sig_sum, axis = 1, thresh = t, methods = 'position')
cands1 = final1.loc[(ss_motif >= 1) & (ss_position >= 1)]
pairs1 = set([(int(x.GRE_ID), x.TF_ID) for _,x in cands1.iterrows()])
ss_motif = final2.apply(sig_sum, axis = 1, thresh = t, methods = 'motif')
ss_position = final2.apply(sig_sum, axis = 1, thresh = t, methods = 'position')
cands2 = final2.loc[(ss_motif >= 1) & (ss_position >= 1)]
pairs2 = set([(int(x.GRE_ID), x.TF_ID) for _,x in cands2.iterrows()])
overlap = pairs1.intersection(pairs2)
venn.venn2([pairs1,pairs2])
diffy = pairs2.difference(pairs1)

dbo = different_best_ones
dbo2 = different_best_ones2
for pair in random.sample(diffy,10):
    GRE_ID = pair[0]
    TF_ID = pair[1]
    e1 = dbo.loc[(dbo.GRE_ID == GRE_ID)&(dbo.TF_ID == TF_ID)]
    e2 = dbo2.loc[(dbo2.GRE_ID == GRE_ID)&(dbo2.TF_ID == TF_ID)]
    ps1 = e1.peaks.iloc[0]
    Tms1 = e1.TF_shape.iloc[0]
    TmI1 = int(e1.TF_motif_ID.iloc[0])
    ps2 = e2.peaks.iloc[0]
    Tms2 = e2.TF_shape.iloc[0]
    TmI2 = int(e2.TF_motif_ID.iloc[0])
    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF_ID, peak_set = ps1, TF_motif_shape = Tms1, TF_motif_ID = TmI1)
    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF_ID, peak_set = ps2, TF_motif_shape = Tms2, TF_motif_ID = TmI2)
 
#decision: choice of TFmotif-GRE match based on all 3 motif pvalues
    
final = final2
    
#V.##---possible criterion: how many methods show significance based on thresholds?-------------####

#standard alpha threshold
alpha = 0.05
#Bonferroni corrected threshold
bf = alpha/(len(all_GREs)*len(all_TFs))
#threshold based on 1 expected significant hit by chance
#with n = number of attempts to find a match = len(all_GREs)*len(all_TFs)
#and p = significance threshold for mapping = probability to find a same/higher scoring match by chance = p_motif*p_position assuming that scores for both methods are independent
#If we choose one threshold for both methods:
# E = p^2*n -> p = sqrt(E/n)
#for one expected hit by chance, we define the threshold as  
pE1 = np.sqrt(1/len(final))    

#try to define another threshold by forcing every TF to have at least one match

#see XII for analysis of coverage and connectivity based on pvalue thresholds

final['ss_pos_alpha'] = final.apply(sig_sum, axis = 1, thresh = alpha, methods = "position")
final['ss_motif_alpha'] = final.apply(sig_sum, axis = 1, thresh = alpha, methods = "motif")
final['ss_pos_pE1'] = final.apply(sig_sum, axis = 1, thresh = pE1, methods = "position")
final['ss_motif_pE1'] = final.apply(sig_sum, axis = 1, thresh = pE1, methods = "motif")
final['ss_pos_bf'] = final.apply(sig_sum, axis = 1, thresh = bf, methods = "position")
final['ss_motif_bf'] = final.apply(sig_sum, axis = 1, thresh = bf, methods = "motif")
final['ss_alpha'] = final.ss_pos_alpha + final.ss_motif_alpha
final['ss_pE1'] = final.ss_pos_pE1 + final.ss_motif_pE1
final['ss_bf'] = final.ss_pos_bf + final.ss_motif_bf


#filtering for hiqh confidence matches based on above thresholds

pE1_filtered = final.loc[(final.ss_pos_pE1>= 1) & (final.ss_motif_pE1>= 1)]
alpha_filtered = final.loc[(final.ss_pos_alpha>= 1) & (final.ss_motif_alpha>= 1)]
bf_filtered = final.loc[(final.ss_pos_bf>= 1) & (final.ss_motif_bf>= 1)]

#VI.##---------------------------correlations between different mapping scores------------------####

variants = ['raw', 'unique_TF', 'p_alpha_filtered', 'p_pE1_filtered', 'p_bf_filtered']
corrs = {k:{} for k in variants}
for k,d in zip(variants, [complete_result, final, alpha_filtered, pE1_filtered, bf_filtered]):
    corrs[k]['all_cols'] = d.corr(method = 'spearman')
    corrs[k]['pval_cols'] = d[['pv', 'pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].corr(method = 'spearman')
c = corrs["raw"]['pval_cols']
c2 = corrs['unique_TF']['pval_cols']
c3 = corrs['p_alpha_filtered']['pval_cols']
c4 = corrs['p_pE1_filtered']['pval_cols']
c5 = corrs['p_bf_filtered']['pval_cols']
c.columns = c2.columns = c3.columns = c4.columns = c5.columns = ['tomtom uncut','tomtom GRE cut','tomtom both cut','position genes','position hits','position nucleotides', "position peak center"]
c.index = c2.index = c3.index = c4.index = c5.index = c.columns
with PdfPages('correlations_methods.pdf') as pp:
    fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
    ax = fig.add_subplot(221)
    sns.heatmap(c2,cmap = "Blues", annot = True, vmin = 0, vmax = 1)
    ax.set_title("correlation of mapping pvalues \n using all combinations \n $n_{maps}$ = %d"%len(final) )
    ax = fig.add_subplot(222)
    sns.heatmap(c3,cmap = "Blues", annot = True, vmin = 0, vmax = 1)
    ax.set_title("correlation of mapping pvalues \n using only combinations with $pv_{all}$ < 0.05  \n $n_{maps}$ = %d"%len(alpha_filtered))
    ax = fig.add_subplot(223)
    sns.heatmap(c4,cmap = "Blues", annot = True, vmin = 0, vmax = 1)
    ax.set_title("correlation of mapping pvalues \n using only combinations with $pv_{all}$ < %1.1e  \n $n_{maps}$ = %d"%(pE1, len(pE1_filtered)))
    ax = fig.add_subplot(224)
    sns.heatmap(c5,cmap = "Blues", annot = True, vmin = 0, vmax = 1)
    ax.set_title("correlation of mapping pvalues \n using only combinations with $pv_{all}$ < %1.1e  \n $n_{maps}$ = %d"%(bf, len(bf_filtered)))
    plt.tight_layout()
    pp.savefig()



        
#VII.#--------number of methods  indicating significance based on different thresholds----------####

from brokenaxes import brokenaxes
with PdfPages('no_sig_methods.pdf') as pp:
    fig = plt.figure(figsize=(16, 8.27), dpi=300)
    sps1, sps2, sps3 = mpl.gridspec.GridSpec(1,3)
    #alpha
    bax = fig.add_subplot(131)
    bax.bar(x = [0,1,2,3,4,5,6,7], height = final.groupby(by = 'ss_alpha').count().GRE_ID, color = ['#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'])
    bax.set_title('Number of matches supported by n methods \n based on pvalue threshold {}'.format(alpha))
    bax.set_ylabel('# matches')
    bax.yaxis.set_label_coords(-0.15, 0.3)
    bax.set_xlabel('# methods supporting match')
    bax.text(x = 5.2, y = 400, s = "n >= 5: {:d}".format(len(final.loc[final.ss_alpha >= 5])))
    #pE1
    bax = brokenaxes(ylims = ((0,4000),(12500,13000)), subplot_spec = sps2, wspace = 0.2, despine = False, d = 0.01)
    bax.bar(x = [0,1,2,3,4,5,6,7], height = final.groupby(by = 'ss_pE1').count().GRE_ID, color = ['#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'])
    bax.set_title('Number of matches supported by n methods \n based on pvalue threshold {:.1e}'.format(pE1))
    bax.set_ylabel('# matches')
    bax.big_ax.yaxis.set_label_coords(-0.15, 0.3)
    bax.set_xlabel('# methods supporting match')
    bax.text(x = 5.2, y = 250, s = "n >= 5: {:d}".format(len(final.loc[final.ss_pE1 >= 5])))
    #bf
    bax = brokenaxes(ylims = ((0,100), (500,700),(18200, 18500)), subplot_spec = sps3, wspace = 0.2, despine = False, d = 0.01)
    bax.bar(x = [0,1,2,3,4,5,6,7], height = final.groupby(by = 'ss_bf').count().GRE_ID, color = ['#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'])
    bax.set_title('Number of matches supported by n methods \n based on pvalue threshold {:.1e}'.format(bf))
    bax.set_ylabel('# matches')
    bax.big_ax.yaxis.set_label_coords(-0.15, 0.3)
    bax.set_xlabel('# methods supporting match')
    bax.text(x = 5.2, y = 25, s = "n >= 5: {:d}".format(len(final.loc[final.ss_bf >= 5])))
    #plt.tight_layout()
    pp.savefig()
    

  
#VIII.---dictionaries containing mapping subsets based on different thresholds for each method--####

methods = ['uncut', 'GRE_cut', 'both_cut', 'genes', 'hits', 'nucleotides', 'peak_center']
filtered_alpha = {}
for m in methods:
    f = filter_pvals(final, methods = m, thresh = alpha)
    filtered_alpha[m] = set([(int(x.GRE_ID), x.TF_ID) for _, x in f.iterrows()])
    
filtered_bf = {}
for m in methods:
    f = filter_pvals(final, methods = m, thresh = bf)
    filtered_bf[m] = set([(int(x.GRE_ID), x.TF_ID) for _, x in f.iterrows()])
    
filtered_pE1 = {}
for m in methods:
    f = filter_pvals(final, methods = m, thresh = pE1)
    filtered_pE1[m] = set([(int(x.GRE_ID), x.TF_ID) for _, x in f.iterrows()])



#IX.-------------------sets and subsets of significnace for different methods-------------------####

intersec_motif = filtered_pE1['uncut'].intersection(filtered_pE1['GRE_cut']).intersection(filtered_pE1['both_cut'])
union_motif = filtered_pE1['uncut'].union(filtered_pE1['GRE_cut']).union(filtered_pE1['both_cut'])
intersec_position = filtered_pE1['genes'].intersection(filtered_pE1['hits']).intersection(filtered_pE1['nucleotides'])
union_position = filtered_pE1['genes'].union(filtered_pE1['hits']).union(filtered_pE1['nucleotides'])

import matplotlib_venn as venn
with PdfPages('venns_confidence_sets_pE1.pdf') as pp:
    fig = plt.figure(figsize = (11.69, 8.27), dpi=300)
    fig.suptitle("based on threshold p <= {:1.1e}".format(pE1))
    ax = fig.add_subplot(221)
    #the sets based on motif_mapping
    venn.venn3([filtered_pE1['uncut'], filtered_pE1['GRE_cut'], filtered_pE1['both_cut']], set_labels= ['uncut', 'GRE_cut', 'both_cut'])
    ax.set_title('significant motif-based matches')
    ax = fig.add_subplot(222)
    #the sets based on gene/hit/nucleotide mapping
    venn.venn3([filtered_pE1['genes'], filtered_pE1['hits'], filtered_pE1['nucleotides']], set_labels= ['genes', 'hits', 'nucleotides'])
    ax.set_title('significant position-based matches')
    ax = plt.subplot2grid((2,2), (1,0), colspan = 2)
    venn.venn2([intersec_motif, intersec_position], set_labels = ['intersection motif','intersection position'])
    ax.set_title('overlap of motif- and position-based matches')
    pp.savefig()


#grande_union = union_position.union(union_motif)
#grande_intersec = intersec_motif.intersection(intersec_position)


#TF-GRE pairs based on different numbers of significance methods
quality_sets_bf = {}
quality_sets_pE1 = {}
for i in [0,1,2,3,4,5,6,7]:
    quality_sets_bf[i] = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[final.ss_bf == i].iterrows()])
    quality_sets_pE1[i] = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[final.ss_pE1 == i].iterrows()])

#matches that were not significant before cutting but significant after
signified_set_bf = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[(final.pv >= bf) & ((final.pv_cut_GRE <= bf) | (final.pv_cut_both <= bf))].iterrows()])
signified_set_pE1 = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[(final.pv >= pE1) & ((final.pv_cut_GRE <= pE1) | (final.pv_cut_both <= pE1))].iterrows()])

#reverse: significant before, not after
designified_set_bf = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[(final.pv <= bf) & ((final.pv_cut_GRE >= bf) | (final.pv_cut_both >= bf))].iterrows()])
designified_set_pE1 = set([(int(x.GRE_ID), x.TF_ID) for _,x in final.loc[(final.pv <= pE1) & ((final.pv_cut_GRE >= pE1) | (final.pv_cut_both >= pE1))].iterrows()])


#visualize motif mappings for different sets
from tomtom_for_candidates import tomtom_candidates


#this one works with mapping data with unique TF-GRE matches
for pair in quality_sets_pE1[4]:
    e = final2.loc[(final2.GRE_ID == pair[0])&(final2.TF_ID == pair[1])]
    GRE_ID = e.GRE_ID.iloc[0]
    TF_ID = e.TF_ID.iloc[0]
    peak_set = e.peaks.iloc[0]
    TF_motif_shape = e.TF_shape.iloc[0]
    TF_motif_ID = e.TF_motif_ID.iloc[0]
    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF_ID, peak_set = peak_set, TF_motif_shape = TF_motif_shape, TF_motif_ID = TF_motif_ID)

#this one works with mapping data with all TF-GRE matches 
for pair in njap:
    e = complete_result.loc[(df.GRE_ID == pair[0])&(df.TF_ID == pair[1])]
    GRE_IDs = e.GRE_ID
    TF_IDs = e.TF_ID
    peak_sets = e.peaks
    TF_motif_shapes = e.TF_shape
    TF_motif_IDs = e.TF_motif_ID
    for GRE_ID, TF_ID, peak_set, TF_motif_shape, TF_motif_ID in zip(GRE_IDs, TF_IDs, peak_sets, TF_motif_shapes, TF_motif_IDs):
        tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF_ID, peak_set = peak_set, TF_motif_shape = TF_motif_shape, TF_motif_ID = TF_motif_ID)




#X.------------------------------analysis of position-based mapping-------------------------#######

df = ml_positions


#plot pvalues of all position-based matches
pp = PdfPages('pval_genes_vs_pval_hits_position_mapping.pdf')
alpha = 0.05
bf = alpha/len(df) 
all_hits = df.loc[(df.pv_genes < bf) & (df.pv_hits < bf) & (df.pv_nucleotides < bf)]
fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
ax = fig.add_subplot(2,2,1)
ax.plot(df.pv_genes, df.pv_hits, '.',markersize = 1, label = "p values",  color = '#9ecae1')
ax.plot(all_hits.pv_genes, all_hits.pv_hits, '.', markersize = 2, label = 'significant in all mappings', color = '#3182bd')
ax.set_xlabel('pvalue based on genes')
ax.set_ylabel('pvalue based on hits')
ax.vlines(alpha, ymin = min(df.pv_hits.loc[df.pv_hits != min(df.pv_hits)]), ymax = 1.5, colors = "r", linestyles = 'dashed', linewidth = 1, label = '$p = \\alpha = %s$'%alpha)
ax.vlines(bf, ymin = min(df.pv_hits.loc[df.pv_hits != min(df.pv_hits)]), ymax = 1.5, colors = "r", linestyles = 'dotted', linewidth = 1, label = '$p = \\alpha_{BF}$ = %3.2e'%bf)
ax.hlines([alpha,bf], xmin = min(df.pv_genes), xmax = 1.5, colors = "r", linestyles = ['dashed', 'dotted'], linewidth = [1,1])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([1e-40, 5])
ax.set_ylim([1e-40, 5])
n_alpha_genes = len(df.loc[df.pv_genes < alpha])
n_alpha_positions = len(df.loc[df.pv_hits < alpha])
n_alpha_both = len(df.loc[(df.pv_genes < alpha) & (df.pv_hits < alpha)])
n_bf_genes = len(df.loc[df.pv_genes < bf])
n_bf_positions = len(df.loc[df.pv_hits < bf])
n_bf_both = len(df.loc[(df.pv_genes < bf) & (df.pv_hits < bf)])
ax.text(x = 0.1, y = 0.5, s = ' {} matches below gene threshold, \n {} matches below hits threshold, \n {} matches below both thresholds'.format(n_bf_genes, n_bf_positions, n_bf_both), transform = ax.transAxes)
ax.set_title('pvalues for matches based on gene mapping and hit mapping')
ax.legend(loc = 3)

ax = fig.add_subplot(2,2,2)
ax.plot(df.pv_genes, df.pv_nucleotides, '.',markersize = 1, label = "p values", color = '#9ecae1')
ax.plot(all_hits.pv_genes, all_hits.pv_nucleotides, '.', markersize = 2, label = 'significant in all mappings', color = '#3182bd')
ax.set_xlabel('pvalue based on genes')
ax.set_ylabel('pvalue based on nucleotide hits')
ax.vlines(alpha, ymin = min(df.pv_nucleotides.loc[df.pv_nucleotides != min(df.pv_nucleotides)]), ymax = 1.5, colors = "r", linestyles = 'dashed', linewidth = 1, label = '$p = \\alpha = %s$'%alpha)
ax.vlines(bf, ymin = min(df.pv_nucleotides.loc[df.pv_nucleotides != min(df.pv_nucleotides)]), ymax = 1.5, colors = "r", linestyles = 'dotted', linewidth = 1, label = '$p = \\alpha_{BF}$ = %3.2e'%bf)
ax.hlines([alpha,bf], xmin = min(df.pv_genes), xmax = 1.5, colors = "r", linestyles = ['dashed', 'dotted'], linewidth = [1,1])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([1e-40, 5])
ax.set_ylim([1e-40, 5])
n_alpha_genes = len(df.loc[df.pv_genes < alpha])
n_alpha_positions = len(df.loc[df.pv_nucleotides < alpha])
n_alpha_both = len(df.loc[(df.pv_genes < alpha) & (df.pv_nucleotides < alpha)])
n_bf_genes = len(df.loc[df.pv_genes < bf])
n_bf_positions = len(df.loc[df.pv_nucleotides < bf])
n_bf_both = len(df.loc[(df.pv_genes < bf) & (df.pv_nucleotides < bf)])
ax.text(x = 0.1, y = 0.5, s = ' {} matches below gene threshold, \n {} matches below nucleotide threshold, \n {} matches below both thresholds'.format(n_bf_genes, n_bf_positions, n_bf_both), transform = ax.transAxes)
ax.set_title('pvalues for matches based on gene mapping and nucleotide mapping')
ax.legend(loc = 3)

ax = fig.add_subplot(2,2,3)
ax.plot(df.pv_hits, df.pv_nucleotides, '.',markersize = 1, label = "p values", color = '#9ecae1')
ax.plot(all_hits.pv_hits, all_hits.pv_nucleotides, '.', markersize = 2, label = 'significant in all mappings', color = '#3182bd')
ax.set_xlabel('pvalue based on hits')
ax.set_ylabel('pvalue based on nucleotide hits')
ax.vlines(alpha, ymin = min(df.pv_nucleotides.loc[df.pv_nucleotides != min(df.pv_nucleotides)]), ymax = 1.5, colors = "r", linestyles = 'dashed', linewidth = 1, label = '$p = \\alpha = %s$'%alpha)
ax.vlines(bf, ymin = min(df.pv_nucleotides.loc[df.pv_nucleotides != min(df.pv_nucleotides)]), ymax = 1.5, colors = "r", linestyles = 'dotted', linewidth = 1, label = '$p = \\alpha_{BF}$ = %3.2e'%bf)
ax.hlines([alpha,bf], xmin = min(df.pv_hits.loc[df.pv_hits != min(df.pv_hits)]), xmax = 1.5, colors = "r", linestyles = ['dashed', 'dotted'], linewidth = [1,1])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([1e-40, 5])
ax.set_ylim([1e-40, 5])
n_alpha_hits = len(df.loc[df.pv_hits < alpha])
n_alpha_positions = len(df.loc[df.pv_nucleotides < alpha])
n_alpha_both = len(df.loc[(df.pv_hits < alpha) & (df.pv_nucleotides < alpha)])
n_bf_hits = len(df.loc[df.pv_hits < bf])
n_bf_positions = len(df.loc[df.pv_nucleotides < bf])
n_bf_both = len(df.loc[(df.pv_hits < bf) & (df.pv_nucleotides < bf)])
ax.text(x = 0.1, y = 0.5, s = ' {} matches below hits threshold, \n {} matches below nucleotide threshold, \n {} matches below both thresholds'.format(n_bf_hits, n_bf_positions, n_bf_both), transform = ax.transAxes)
ax.set_title('pvalues for matches based on hits mapping and nucleotide mapping')
ax.legend(loc = 3)

ax = fig.add_subplot(2,2,4)
n_all_hits = len(df.loc[(df.pv_genes < bf) & (df.pv_hits < bf) & (df.pv_nucleotides < bf)])
ax.text(x = 0.5, y = 0.5, s = 'matches significant based on all three mappings: {}'.format(n_all_hits), transform = ax.transAxes, horizontalalignment='center', verticalalignment='center')
ax.set_axis_off()

plt.tight_layout()
pp.savefig()
pp.close()


#XI.--------store mapping list as csv with only columns used for evaluation------------------#######

final_cur = final[['GRE_ID', 'TF_ID', 'pv', 'pv_cut_GRE', 'pv_cut_both', 'pv_genes',
                   'pv_hits', 'pv_nucleotides', 'pv_peak_center', 'ss_pos_pE1', 'ss_motif_pE1',
                   'ss_pos_bf', 'ss_motif_bf', 'ss_pos_alpha', 'ss_motif_alpha', 'ss_alpha',
                   'ss_pE1', 'ss_bf']]

final.to_csv(open('final_mapping_result.csv', 'w'))
final_cur.to_csv(open('cured_final_mapping_result.csv', 'w'))



#XII.------iterate over multiple thresholds for influence on coverage/connectivity ----------#######


county = pb.ProgressBar()
connectivity = {}
coverage = {}
no_connections = {}
all_TFs = set(final.TF_ID)
all_GREs = set(final.GRE_ID)
for t in county([1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 
                 0.4, 0.3, 0.2, 0.1, 0.09, 0.09, 0.07, 0.06, 0.05, 
                 0.04, 0.03, 0.02, 0.01, pE1, 0.005, 0.001, 0.0005, 
                 0.0001, 1e-05,1e-06, 1e-07, 1e-08, 1e-09, 1e-10, 
                 1e-11, 1e-12]):#, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17]):
    ss_motif = final.apply(sig_sum, axis = 1, thresh = t, methods = 'motif')
    #ss_position = final.apply(sig_sum, axis = 1, thresh = t, methods = 'position')
    #ss_both = ss_motif + ss_position
    #cands = final.loc[(ss_motif >= 3) & (ss_position >= 2)]# & (ss_both >=5)]
    cands = final.loc[(ss_motif >= 3)]# & (ss_position >= 2)]# & (ss_both >=5)]
    no_connections[t] = len(cands)
    TFs = set(cands.TF_ID)
    GREs = set(cands.GRE_ID)
    connectivity[t] = {}
    connectivity[t]['TFs'] = len(cands)/len(TFs)
    connectivity[t]['GREs'] = len(cands)/len(GREs)
    connectivity[t]['both'] = len(cands)*2/(len(TFs)+(len(GREs)))
    coverage[t] = {}
    coverage[t]['TFs'] = len(TFs)/len(all_TFs)
    coverage[t]['GREs'] = len(GREs)/len(all_GREs)
    
lists = sorted(coverage.items())
ps, coverages = zip(*lists)
lists = sorted(connectivity.items())
ps, connectivities = zip(*lists)
lists = sorted(no_connections.items())
ps, no_connections = zip(*lists)
#
#
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}
axe = {'titlesize' : 24, 
       'labelsize' : 18}
       
line = {'linewidth' : 3,
        'markersize' : 10}

plt.rc('axes', **axe)
plt.rc('lines', **line)
plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)
plt.rc('legend', fontsize = 'x-large')
#plt.style.use('default')
#plt.style.use('presentation')
#plt.rc('font', **font)
#with PdfPages('motif_3_pval_vs_coverage_and_connectivity.pdf') as pp:
fig = plt.figure(figsize = (16, 12), dpi = 300)
ax = fig.add_subplot(311)
ax.plot(list(ps),list(no_connections), '-', color = "k",linestyle = 'dashed', label = "#")
ax.set_xscale('log')
ax.set_yscale("log")
#ax.set_xlabel('motif p-value threshold')
ax.set_xticklabels([])
ax.set_ylabel('# total connections')
ax.vlines(alpha, ymin = 0, ymax = 1, colors = '#fc8d59', label = r'$\alpha$ = 0.05')
#ax.vlines(pE1, ymin = 0, ymax = 1, colors = '#e34a33', linestyles = 'dashed', label = r'$pE1$ = %1.1e'%pE1)
ax.vlines(bf, ymin = 0, ymax = 1, colors = '#b30000',  label = r'$\alpha_{BF}$ = %1.1e'%bf)
ax.set_title('Total number')
ax.grid()
ax = fig.add_subplot(312)
ax.plot(list(ps),[x['TFs'] for x in coverages], '-', linestyle = 'dashed',color = "#7fc97f", label = "TFs")
ax.plot(list(ps),[x['GREs'] for x in coverages], '-',linestyle = 'dashed', color = "#7b3294", label = 'GREs')
ax.set_xscale('log')
#ax.set_xlabel('motif p-value threshold')
ax.set_xticklabels([])
ax.set_ylabel('ratio of mapped TFs and GREs')
#    ax.hlines(126, xmin =1e-14, xmax = 1, colors = "#7fc97f", linestyles = "dashed", label = '# TFs total' )
#    ax.hlines(152, xmin =1e-14, xmax = 1, colors = "#7b3294", linestyles = "dashed", label = '# GREs total')
ax.vlines(alpha, ymin = 0, ymax = 1, colors = '#fc8d59', label = r'$\alpha$ = 0.05')
#ax.vlines(pE1, ymin = 0, ymax = 1, colors = '#e34a33', linestyles = 'dashed', label = r'$pE1$ = %1.1e'%pE1)
ax.vlines(bf, ymin = 0, ymax = 1, colors = '#b30000', label = r'$\alpha_{BF}$ = %1.1e'%bf)
ax.set_title('Coverage')
ax.grid()
ax.legend(loc = 2)
ax = fig.add_subplot(313)
ax.plot(list(ps), [x["TFs"] for x in connectivities], '-', color = "#7fc97f",linestyle = 'dashed', label = "GREs/TF")
ax.plot(list(ps), [x["GREs"] for x in connectivities], '-', color = "#7b3294",linestyle = 'dashed', label = "TFs/GRE")
ax.plot(list(ps), [x["both"] for x in connectivities], '-', color = "k",linestyle = 'dashed', label = "combined")
ax.set_xscale('log')
ax.set_xlabel('p-value threshold')
ax.set_ylabel('average connections per item')
ax.set_ylim([0,10])
ax.vlines(alpha, ymin = 0, ymax = 152, colors = '#fc8d59', label = r'$\alpha$ = 0.05')
#ax.vlines(pE1, ymin = 0, ymax = 152, colors = '#e34a33', linestyles = 'dashed', label = r'$pE1$ = %1.1e'%pE1)
ax.vlines(bf, ymin = 0, ymax = 152, colors = '#b30000', label = r'$\alpha_{BF}$ = %1.1e'%bf)
ax.set_title('TFs/GRE and GREs/TFs')
ax.legend(loc = 2)
ax.grid()
plt.tight_layout()
plt.savefig('motif_3_pval_vs_no_and_coverage_and_connectivity.eps')
#
#
#
#
#plt.hist(np.log10(hits.evalue_GRE), bins = 100)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.hist(np.log10(hits.evalue_TF.loc[hits.evalue_TF != 0]), bins = 100)
#
#
#
##XIII.-------------------------extract hits based on conditions described below----------------#######
#
## I. Every method has a significance pvalue threshold of 2.6e-6, a value ~ corresponding to one expected false positive match given random background.
#
##II. A match is considered significant if at least one motif-based and one position-based method shows significance AND at least 4/7 total methods are significant.
#
#final['hit'] = (final.ss_motif_pE1 >= 1) & (final.ss_pos_pE1 >= 1) & (final.ss_pE1 >= 4)
#hits = final.loc[final.hit]
#
#hit_GREs = set(hits.GRE_ID)
#hit_TFs = set(hits.TF_ID)
#connectivity = {'TFs': len(hits)/len(hit_TFs), 'GREs':len(hits)/len(hit_GREs)}
#coverage = {'TFs': len(hit_TFs)/len(all_TFs), 'GREs': len(hit_GREs)/len(all_GREs)}
#connections_per_TF = hits.TF_ID.value_counts().to_dict()
#connections_per_GRE = hits.GRE_ID.value_counts().to_dict()
#
#
##new entry for every match: how many TF hits has the involved GRE, how many GRE hits the involved TF?
#hits['no_TF_connections'] = [nan]*len(hits)
#hits['no_GRE_connections'] = [nan]*len(hits) 
#for i,x in hits.iterrows():
#    k = x.TF_ID
#    hits.at[i, 'no_TF_connections'] = int(connections_per_TF[k])
#    k = x.GRE_ID
#    hits.at[i, 'no_GRE_connections'] = int(connections_per_GRE[k])
#    
##correlations(for analysis of number of connections)
#c = hits.corr('spearman')
#
#
#
#
#
##little excursion: Rv1033c binds 14! GREs significantly based on conditions described in XIII
#hits_1033c = hits.loc[hits.TF_ID == 'Rv1033c']
#Rv1033_pairs = set([(int(x.GRE_ID), x.TF_ID) for _,x in hits_1033c.iterrows()])
#
##check motif matching
#for pair in Rv1033_pairs:
#    e = hits.loc[(hits.GRE_ID == pair[0])&(hits.TF_ID == pair[1])]
#    GRE_ID = e.GRE_ID.iloc[0]
#    TF_ID = e.TF_ID.iloc[0]
#    peak_set = e.peaks.iloc[0]
#    TF_motif_shape = e.TF_shape.iloc[0]
#    TF_motif_ID = e.TF_motif_ID.iloc[0]
#    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF_ID, peak_set = peak_set, TF_motif_shape = TF_motif_shape, TF_motif_ID = TF_motif_ID)
#    
#
#
#
#
##XIV.--overlap of TFbound/diffex genes and GRE-associated genes for GREs associated with TFs--#######
#
#
##GRE-gene associations from csv files ('gene_all_gres')
#dir_ = '/Users/rkoch/Documents/Data_and_Scripts/Data/gene_all_gres/'
#GRE_genes = {}
#for GRE in all_GREs:
#    GRE_genes[GRE] = set()
#    ID = GRE+1
#    file = dir_ + 'GRE{}.csv'.format(ID)
#    if os.path.exists(file):
#        with open(file, 'r') as f:
#            lines = [x.split('"') for x in f.readlines()]
#            genes = [x[1] for x in lines[1:] if len(x)>1]
#            GRE_genes[GRE] = set(genes)
#
##GRE-gene associations from GRE-corem--corem-genes associations
#GRE_cor_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/GRE_cor_genes.p', 'rb'))
#
##TF-gene associations (TFOE)
#TF_diffex_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts//Dictionaries/TF_diffex_dic.p', 'rb'))
#counywouny = pb.ProgressBar()
#
##here we find for each TF all GREs mapped to it, sort those GREs according to the median pv of the 7 methods, 
##then stepwise include more and more GREs and associate the TF with the genes associated to the GREs. Then, 
##recall of genes diffex upon TFOE via the set of TF-GRE-(corem) associated genes is evaluated as 
##hypergeom. pvalue of the overlap of both gene sets 
#TF_recall_diffex_with_GREs = {}
#for TF in counywouny(all_TFs):
#    TF_recall_diffex_with_GREs[TF] = []
#    if TF in TF_diffex_dic:
#        TFOE_genes = TF_diffex_dic[TF]
#    else:
#        TFOE_genes = []
#    n_TFOE_genes = len(TFOE_genes)
#    df = final.loc[final.TF_ID == TF]
#    df['meds'] = df[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
#    df.sort_values(by = 'meds', inplace = True)
#    TF_assoc_genes = set()      
#    for i,row in zip(range(200), df.iterrows()):
#        med_pv = row[1].meds
#        if med_pv > 0.05:
#            break
#        GRE = row[1].GRE_ID
#        if GRE in GRE_genes:                        #or GRE in GRE_genes:
#            TF_assoc_genes.update(GRE_genes[GRE])   #or .update(GRE_genes[GRE])
#            n_TF_assoc_genes = len(TF_assoc_genes)
#            n_both = len(TF_assoc_genes.intersection(TFOE_genes))
#            N = 3996 #number of MTB genes
#            hpd = ss.hypergeom(N, n_TFOE_genes, n_TF_assoc_genes) #hypergeometric distribution of drawing n_TF_assoc_genes from N total genes of which n_TFOE_genes are class I
#            pv = hpd.sf(n_both-1) #p(n >= n_both); is identical to 1-cdf(kboth-1); sf = survival function = 1-cdf
#            TF_recall_diffex_with_GREs[TF].append([i, med_pv, pv,])
#        
#        
#with PdfPages('overlap_TF_diffex_genes_with_TF_GRE_genes_based_on_no_of_mapped_GREs.pdf') as pp:
#    fig = plt.figure(figsize = (116.9, 82.7))
#    i=1
#    for TF, curve in TF_recall_diffex_with_GREs.items():
#        if curve:
#            ax = fig.add_subplot(11,11,i)
#            df = pd.DataFrame(curve, columns = ['index', 'median_pv', 'overlap_pv'])
#            ax.plot(df.index, df.overlap_pv, 'x-')
#            ax.set_ylim(-0.1,1.1)
#            #ax.set_yscale("log")
#            ax.set_xlim(0,10)
#            ax.set_title('{}'.format(TF))
#            i+=1
#    pp.savefig()
#
#  
##The same procedure with set of TF-bound genes (instead of diffex genes)
#
#TF_bound_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts//Dictionaries/pub_chip_dic.p', 'rb'))
#
#counywouny = pb.ProgressBar()
#TF_recall_bound_with_GREs = {}
#for TF in counywouny(all_TFs):
#    TF_recall_bound_with_GREs[TF] = []
#    if TF in TF_bound_dic:
#        TF_bound_genes = TF_bound_dic[TF]
#    else:
#        TF_bound_genes = []
#    n_TF_bound_genes = len(TF_bound_genes)
#    df = final.loc[final.TF_ID == TF]
#    df['meds'] = df[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
#    df.sort_values(by = 'meds', inplace = True)
#    TF_assoc_genes = set()      
#    for i,row in zip(range(200), df.iterrows()):
#        med_pv = row[1].meds
#        if med_pv > 0.05:
#            break
#        GRE = row[1].GRE_ID
#        if GRE in GRE_cor_genes:                    #or GRE in GRE_genes:
#            TF_assoc_genes.update(GRE_cor_genes[GRE])   #or .update(GRE_genes[GRE])
#            n_TF_assoc_genes = len(TF_assoc_genes)
#            n_both = len(TF_assoc_genes.intersection(TF_bound_genes))
#            N = 3996 #number of MTB genes
#            hpd = ss.hypergeom(N, n_TF_bound_genes, n_TF_assoc_genes) #hypergeometric distribution of drawing n_TF_assoc_genes from N total genes of which n_TF_bound_genes are class I
#            pv = hpd.sf(n_both-1) #p(n >= n_both); is identical to 1-cdf(kboth-1); sf = survival function = 1-cdf
#            TF_recall_bound_with_GREs[TF].append([i, med_pv, pv,])   
#    
#with PdfPages('overlap_TF_bound_genes_with_TF_GRE_cor_genes_based_on_no_of_mapped_GREs.pdf') as pp:
#    fig = plt.figure(figsize = (116.9, 82.7))
#    i=1
#    for TF, curve in TF_recall_bound_with_GREs.items():
#        if curve:
#            ax = fig.add_subplot(11,11,i)
#            df = pd.DataFrame(curve, columns = ['index', 'median_pv', 'overlap_pv'])
#            ax.plot(df.index, df.overlap_pv, 'x-')
#            ax.set_ylim(-0.1,1.1)
#            #ax.set_yscale("log")
#            ax.set_xlim(0,10)
#            ax.set_title('{}'.format(TF))
#            i+=1
#    pp.savefig()
#
##some TFs have very low overlap pvalues independent of how many GREs are included (e.g. Rv3133c)
## -> because overlap with the first GRE set is very significant and following GREs do not disturb this significance much
#    
#    
#
#
#
#
#####--overlap of TFbound/diffex genes and GRE-associated genes for best GRE associated with TFs--#######
#
#TF_GRE_gene_set_overlaps = {}
#
#counywouny = pb.ProgressBar()
#for TF in counywouny(all_TFs):
#    TF_GRE_gene_set_overlaps[TF] = {}
#    if TF in TF_bound_dic:
#        TF_bound_genes = set(TF_bound_dic[TF])
#    else:
#        TF_bound_genes = set()
#    if TF in TF_diffex_dic:
#        TF_diffex_genes = set(TF_diffex_dic[TF])
#    else:
#        TF_diffex_genes = set()
#    df = final.loc[final.TF_ID == TF]
#    df['meds'] = df[['pv','pv_cut_GRE', 'pv_cut_both', 'pv_genes', 'pv_hits', 'pv_nucleotides', 'pv_peak_center']].apply(np.median, axis = 1)
#    df.sort_values(by = 'meds', inplace = True)
#    GRE = df.iloc[0].GRE_ID
#    if GRE in GRE_cor_genes:
#        TF_GRE_cor_genes = GRE_cor_genes[GRE]
#    else:
#        TF_GRE_cor_genes = set()
#    if GRE in GRE_genes:
#        TF_GRE_genes = GRE_genes[GRE]
#    else:
#        TF_GRE_genes = set()
#    for t in [TF_bound_genes, TF_diffex_genes]:
#        if t == TF_bound_genes:
#            tt = 'bound'
#        else:
#            tt = 'diffex'
#        for g in [TF_GRE_cor_genes, TF_GRE_genes]:
#            if g == TF_GRE_cor_genes:
#                gg = 'cor'
#            else:
#                gg = 'noncor'
#            
#            nt = len(t)
#            ng = len(g)
#            nb = len(t.intersection(g))
#            N = 3996
#            hpd = ss.hypergeom(N, nt, ng) #hypergeometric distribution of drawing n_TF_assoc_genes from N total genes of which n_TF_bound_genes are class I
#            pv = hpd.sf(nb-1)
#            TF_GRE_gene_set_overlaps[TF]['{}_{}'.format(tt, gg)] = [nt,ng,nb,pv]
#            
#
#big_df = pd.DataFrame(columns =['n_TF_genes','n_GRE_genes','n_both_genes','pv', 'TF'] )
#for TF, values in TF_GRE_gene_set_overlaps.items():
#    df = pd.DataFrame.from_dict(values, orient = 'index')
#    df.columns = ['n_TF_genes','n_GRE_genes','n_both_genes','pv']
#    df['TF'] = len(df)*[TF]
#    big_df = big_df.append(df)
#
#big_df.index.name = 'variant'
#big_df = big_df.set_index('TF', append = True)
#
#big_df.n_GRE_genes = big_df.n_GRE_genes.astype(float)
#big_df.n_TF_genes = big_df.n_TF_genes.astype(float)
#big_df.n_both_genes = big_df.n_both_genes.astype(float)
#big_df.sort_values(by = 'pv', ascending = True)
#
#bound_cor_df = big_df.query("variant == 'bound_cor'")
#alphabf = 0.05/126
#overlap_TFs = set(bound_cor_df.query("pv <= {}".format(alphabf)).index.get_level_values('TF').tolist())
#
#best_match_TFs = set(final.loc[final.ss_pE1 == 7].TF_ID.tolist())
#
#
#
#####-----------------------------TF centered filtering------------------------------------------####
#
#
#
#
#
#ms = ['tomtom uncut', 'tomtom GRE cut', 'tomtom both cut',
#      'position genes', 'position hits', 'position nucleotides',
#      'position peak center', 'median all', 'median motif','median position',
#      'best all', 'best motif', 'best position']     
#
#
#recall_TF_bound_genes_with_GREgenes = {}
#recall_TF_diffex_genes_with_GREgenes = {}
#best_TF_maps = {}
#for m in ms:
#    recall_TF_bound_genes_with_GREgenes[m] = {}
#    recall_TF_diffex_genes_with_GREgenes[m] = {}
#    best_TF_maps[m] = filter_best_TF_maps(final, method = m)
#    for i, map_ in best_TF_maps[m].iterrows():
#        TF = map_.TF_ID
#        #TF_bound_genes = set(TF_bound_dic[TF])
#        if TF in TF_bound_dic:
#            TF_bound_genes = set(TF_bound_dic[TF])
#        else:
#            TF_bound_genes = set()
#        if TF in TF_diffex_dic:
#            TF_diffex_genes = set(TF_diffex_dic[TF])
#        else:
#            TF_diffex_genes = set()
#        n_Tb = len(TF_bound_genes)
#        n_Td = len(TF_diffex_genes)
#        
#        GRE = map_.GRE_ID
#        GRE_assoc_genes = GRE_genes[GRE]
#        n_G = len(GRE_assoc_genes)
#        GRE_cor_assoc_genes = GRE_cor_genes[GRE]
#        n_cG = len(GRE_cor_assoc_genes)
#        
#        n_both_b = len(TF_bound_genes.intersection(GRE_assoc_genes))
#        n_both_d = len(TF_diffex_genes.intersection(GRE_assoc_genes))
#        n_both_b_cor = len(TF_bound_genes.intersection(GRE_cor_assoc_genes))
#        n_both_d_cor = len(TF_diffex_genes.intersection(GRE_cor_assoc_genes))
#        
#        N = 3996 #number of MTB genes
#        
#        hpd = ss.hypergeom(N, n_Tb, n_G) 
#        pv1 = hpd.sf(n_both_b-1)
#        
#        hpd = ss.hypergeom(N, n_Tb, n_cG)
#        pv2 = hpd.sf(n_both_b_cor-1)
#        
#        hpd = ss.hypergeom(N, n_Td, n_G)
#        pv3 = hpd.sf(n_both_d-1)
#        
#        hpd = ss.hypergeom(N, n_Td, n_cG)
#        pv4 = hpd.sf(n_both_d_cor-1)
#        
#        recall_TF_bound_genes_with_GREgenes[m][TF] = [pv1,pv2]
#        recall_TF_diffex_genes_with_GREgenes[m][TF] = [pv3,pv4]
# 
#       
#recall_statistics = {}
#for m in ms:
#    recall_statistics[m] = {'ratio_significant_no_cor':nan,'ratio_significant_cor':nan,'significant_TFs_no_cor':[],'significant_TFs_cor':[]}
#    n_cor = 0
#    n_no_cor = 0
#    for TF, pvs in recall_TF_bound_genes_with_GREgenes[m].items():
#    #for TF, pvs in recall_TF_diffex_genes_with_GREgenes[m].items():
#        if pvs[0]<(0.05/126):
#            n_no_cor +=1
#            recall_statistics[m]['significant_TFs_no_cor'].append(TF)
#        if pvs[1]<(0.05/126):
#            n_cor +=1
#            recall_statistics[m]['significant_TFs_cor'].append(TF)
#    recall_statistics[m]['ratio_significant_no_cor'] = n_no_cor/len(all_TFs)
#    recall_statistics[m]['ratio_significant_cor'] = n_cor/len(all_TFs)
#            
##recall_statistics_diffex = recall_statistics.copy()        
#recall_statistics_bound = recall_statistics.copy()        
#    
#for k,v in recall_statistics_bound.items():
#    print(k)
#    print('ratio of TFs with significant overlap to total TFs based on GRE-genes: {}, \n\
# ratio of TFs with significant overlap to total TFs based on GRE-corem-genes: {}'.format(v['ratio_significant_no_cor'], v['ratio_significant_cor']))
#  
##independent of how the best GRE was derived for each TF, only a fraction of 12-18% of the TFs show significant(<0.05) overlap between
##TF-bound genes and genes associated with the best GRE, either through corems or through other means
#    
##similarly, differentially expressed genes overlapped with the GRE-associated genes for only 5-9% of the TFs
#    
#    #-> indication that
#      
#        #either the TF-GRE maps are shitty OR
#        #the GRE-gene maps are shitty OR
#        #we can't recall those gene sets because of the environmental specificity of the gene-TF connections in the experimental setup
#        
#sig_overlapping_TFs_cor = {}
#sig_overlapping_TFs_no_cor = {}
#for m, vs in recall_statistics_bound.items():
#    for TF in vs["significant_TFs_cor"]:
#        if TF in sig_overlapping_TFs_cor:
#            sig_overlapping_TFs_cor[TF] += 1
#        else:
#            sig_overlapping_TFs_cor[TF] = 1
#    for TF in vs["significant_TFs_no_cor"]:
#        if TF in sig_overlapping_TFs_no_cor:
#            sig_overlapping_TFs_no_cor[TF] += 1
#        else:
#            sig_overlapping_TFs_no_cor[TF] = 1
#
#
##compare the sets of TFs where overlap is significant between bound/diffex genes and genes associated with best GRE for that TF 
##to the set of TFs that have high confidence matches based on 
#
#sig_overlap_TF_bound_cor_set = set(sig_overlapping_TFs_cor.keys())
#sig_overlap_TF_bound_nocor_set = set(sig_overlapping_TFs_no_cor.keys())
##sig_overlap_TF_diffex_cor_set = set(sig_overlapping_TFs_cor.keys())
##sig_overlap_TF_diffex_nocor_set = set(sig_overlapping_TFs_no_cor.keys())
#
##45 high quality mappings with 32 unique TFs
#high_conf_TFs = set(hits.loc[hits.ss_pE1 >=6].TF_ID)
#
#venn.venn3([high_conf_TFs, sig_overlap_TF_diffex_cor_set, sig_overlap_TF_bound_cor_set])
##most of the significally overlapping TFs are in the top 32 TFs based on pv filtering
# 
#
#
###filter hits based on top-n approach with subsequent union/intersection
#best_n = {}
#best_n_indices = {}
#for n in [1000,750,500,400,300,200,100,80,60,40,30,20]:
##    best_n[n] = {}
##    for m in ms:
##        best_n[n][m] = filter_best_n(final, n, method = m)
##    best_n_indices[n] = [set(x.index) for x in best_n[n].values()]
#    index_intersec = set.intersection(*best_n_indices[n])
#    index_union = set.union(*best_n_indices[n])
#    best_n_intersec = final.loc[index_intersec]
#    best_n_union = final.loc[index_union]
#    TFs = set(best_n_intersec.TF_ID)
#    GREs = set(best_n_intersec.GRE_ID)
#    connectivity[n] = {}
#    connectivity[n]['TFs'] = len(best_n_intersec)/len(TFs)
#    connectivity[n]['GREs'] = len(best_n_intersec)/len(GREs)
#    connectivity[n]['both'] = len(best_n_intersec)*2/(len(TFs)+(len(GREs)))
#    coverage[n] = {}
#    coverage[n]['TFs'] = len(TFs)/len(all_TFs)
#    coverage[n]['GREs'] = len(GREs)/len(all_GREs)
#    
##bring in plottable format
#lists = sorted(coverage.items())
#ns, coverages = zip(*lists)
#lists = sorted(connectivity.items())
#ns, connectivities = zip(*lists)
#
#
#with PdfPages('top_n_vs_coverage_and_connectivity.pdf') as pp:
#    fig = plt.figure(figsize = (11.69, 8.27), dpi = 300)
#    ax = fig.add_subplot(211)
#    ax.plot(list(ns),[x['TFs'] for x in coverages], '-', color = "#7fc97f", label = "TFs")
#    ax.plot(list(ns),[x['GREs'] for x in coverages], '-', color = "#7b3294", label = 'GREs')
#    #ax.set_xscale('log')
#    #ax.set_xlabel('motif p-value threshold')
#    ax.set_xticks([])
#    ax.set_ylim(0,1.1)
#    ax.set_ylabel('ratio of successfully mapped TFs and GREs')
##    ax.hlines(126, xmin =1e-14, xmax = 1, colors = "#7fc97f", linestyles = "dashed", label = '# TFs total' )
##    ax.hlines(152, xmin =1e-14, xmax = 1, colors = "#7b3294", linestyles = "dashed", label = '# GREs total')
##    ax.vlines(alpha, ymin = 0, ymax = 1, colors = '#fc8d59', linestyles = 'dotted', label = r'$\alpha$ = 0.05')
##    ax.vlines(pE1, ymin = 0, ymax = 1, colors = '#e34a33', linestyles = 'dotted', label = r'$pE1$ = %1.1e'%pE1)
##    ax.vlines(bf, ymin = 0, ymax = 1, colors = '#b30000', linestyles = 'dotted', label = r'$bf$ = %1.1e'%bf)
#    ax.set_title('Coverage for unions of top n')
#    ax.legend(loc = 2)
#    ax = fig.add_subplot(212)
#    ax.plot(list(ns), [x["TFs"] for x in connectivities], '-', color = "#7fc97f", label = "TFs")
#    ax.plot(list(ns), [x["GREs"] for x in connectivities], '-', color = "#7b3294", label = "GREs")
#    ax.plot(list(ns), [x["both"] for x in connectivities], '-', color = "k", label = "combined")
#    #ax.set_xscale('log')
#    ax.set_xlabel('number of top scoring maps included')
#    ax.set_ylabel('connectivity')
#    ax.set_ylim([0,10])
#    #ax.vlines(alpha, ymin = 0, ymax = 152, colors = '#fc8d59', linestyles = 'dotted', label = r'$\alpha$ = 0.05')
#    #ax.vlines(pE1, ymin = 0, ymax = 152, colors = '#e34a33', linestyles = 'dotted', label = r'$pE1$ = %1.1e'%pE1)
#    #ax.vlines(bf, ymin = 0, ymax = 152, colors = '#b30000', linestyles = 'dotted', label = r'$bf$ = %1.1e'%bf)
#    ax.set_title('Connectivity for unions of top n')
#    ax.legend(loc = 2)
#    plt.tight_layout()
#    pp.savefig()
#
#
##check out the pvalue distributions for all methods
#    
#with PdfPages('pval_distributions_lin.pdf') as pp:
#    fig = plt.figure(figsize = (16, 10), dpi=300)
#    fig.suptitle("p-value distributions for different mapping methods")
#    ax = fig.add_subplot(241)
#    ax.hist(complete_result.pv[~np.isnan(complete_result.pv)], bins = 100)
#    ax.set_title('Native tomtom')
#    ax.set_ylabel('#')
#    ax.set_xticks([])
#    ax.set_ylim(1e2,14e3)
#    #ax.set_xlabel('p-value')
#    ax = fig.add_subplot(242)
#    ax.hist(complete_result.pv_cut_GRE[~np.isnan(complete_result.pv_cut_GRE)], bins = 100)
#    ax.set_title('GRE cut tomtom')
#    ax.set_ylim(1e2,14e3)
#    ax.set_yticks([])
#    ax.set_xticks([])
#    #ax.set_ylabel('#')
#    #ax.set_xlabel('p-value')
#    ax = fig.add_subplot(243)
#    ax.hist(complete_result.pv_cut_both[~np.isnan(complete_result.pv_cut_both)], bins = 100)
#    ax.set_title('Both cut tomtom')
#    ax.set_ylim(1e2,14e3)
#    ax.set_yticks([])
#    ax.set_xticks([])
#    #ax.set_ylabel('#')
#    #ax.set_xlabel('p-value')
#    ax = fig.add_subplot(245)
#    ax.hist(complete_result.pv_genes[~np.isnan(complete_result.pv_genes)], bins = 100)
#    ax.set_title('Gene position')
#    ax.set_ylabel('#')
#    ax.set_xlabel('p-value')
#    ax.set_ylim(1e0,3e5)
#    ax = fig.add_subplot(246)
#    ax.hist(complete_result.pv_hits[~np.isnan(complete_result.pv_hits)], bins = 100)
#    ax.set_title('Hits position')
#    #ax.set_ylabel('#')
#    ax.set_ylim(1e0,3e5)
#    ax.set_yticks([])
#    ax.set_xlabel('p-value')
#    ax = fig.add_subplot(247)
#    ax.hist(complete_result.pv_nucleotides[~np.isnan(complete_result.pv_nucleotides)], bins = 100)
#    ax.set_title('Nucleotides position')
#    ax.set_ylim(1e0,3e5)
#    ax.set_yticks([])
#    #ax.set_ylabel('#')
#    ax.set_xlabel('p-value')
#    ax = fig.add_subplot(248)
#    ax.hist(complete_result.pv_peak_center[~np.isnan(complete_result.pv_peak_center)], bins = 100)
#    ax.set_title('Peak center position')
#    ax.set_ylim(1e0,3e5)
#    #ax.set_ylabel('#')
#    ax.set_yticks([])
#    ax.set_xlabel('p-value')
#    #plt.tight_layout()
#    pp.savefig()
#
#
#

####-----------------------------extract regulatory modules-------------------------------------####



#GRE-gene associations from GRE-corem--corem-genes associations
#GRE_cor_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/GRE_cor_genes.p', 'rb'))


#create condition-specific GRE-gene mappings for condition-specific regulatory modules
            
###import static connection data from files as global variables

#GRE-gene connections
#dir_ = '/Users/rkoch/Documents/Data_and_Scripts/Data/gene_all_gres/'
#GRE_genes = {}
#for GRE in all_GREs:
#    GRE_genes[GRE] = set()
#    ID = GRE+1
#    file = dir_ + 'GRE{}.csv'.format(ID)
#    if os.path.exists(file):
#        with open(file, 'r') as f:
#            lines = [x.split('"') for x in f.readlines()]
#            genes = [x[1] for x in lines[1:] if len(x)>1]
#            GRE_genes[GRE] = set(genes)



#GRE-corem connections
GRE_corems = pd.read_csv('/Users/rkoch/Documents/Data_and_Scripts/Data/gre_corem.csv')
GRE_corems.gre = GRE_corems.gre -1#correction for python-realworld index shift
unique_GREs = set(GRE_corems.gre)
GRE_corem_dic = {}
for GRE in unique_GREs:
    corems = list(GRE_corems.loc[GRE_corems.gre == GRE].corem.values)
    GRE_corem_dic[GRE] = corems

#corem-condition connections
corem_cond_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/corem_cond_dic.p', 'rb'))
#corem-gene connections
corems_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/corems_genes.p', 'rb'))

#use generating function to extract GRE-gene connections filtered for all condition
all_conds = set([x for y in corem_cond_dic.values() for x in y])
cond_GRE_genes = {}
for cond in all_conds:
    cond_GRE_genes[cond] = GRE_genes_with_condition(cond)   
cond_GRE_genes['unconditional'] =  GRE_genes_with_condition("all")   


#choose a significant set based on visual inspection of connectivity and coverage
p = 1e-3
ss_motif = final.apply(sig_sum, axis = 1, thresh = p, methods = 'motif')
hits = final.loc[ss_motif == 3]

    #alternative: choose the best GRE match for every TF
hits = pd.DataFrame(columns = final.columns)
for TF in set(final.TF_ID):
    df = final.loc[final.TF_ID == TF]
    df.sort_values('pv', inplace = True)
    hits = hits.append(df.iloc[0])

    #alternative: mega streng, high confidence 

ss_motif = final.apply(sig_sum, axis = 1, thresh = pE1, methods = 'motif')
ss_position = final.apply(sig_sum, axis = 1, thresh = pE1, methods = 'position')
ss_both =  final.apply(sig_sum, axis = 1, thresh = pE1, methods = "both")
strenge_hits = final.loc[(ss_motif==3) & (ss_position >=3)]

#create condition-specific TF-gene mapping via GREs, corems and conditions
all_TFs = set(hits.TF_ID)
cond_TF_GRE_genes = {}
counywouny = pb.ProgressBar()
for TF in counywouny(all_TFs):
    df = hits.loc[hits.TF_ID == TF]
    GREs = df.GRE_ID
    for GRE in GREs:
        if GRE in cond_GRE_genes['hypoxia']:
            if TF in cond_TF_GRE_genes:
                cond_TF_GRE_genes[TF].update(set(cond_GRE_genes['hypoxia'][GRE]))
            else:
                cond_TF_GRE_genes[TF] = set(cond_GRE_genes['hypoxia'][GRE])
                
#create unconditional TF-GRE-gene mapping 
all_TFs = set(hits.TF_ID)
uncond_TF_GRE_genes = {}
counywouny = pb.ProgressBar()
for TF in counywouny(all_TFs):
    df = hits.loc[hits.TF_ID == TF]
    GREs = df.GRE_ID
    for GRE in GREs:
        if GRE in cond_GRE_genes['unconditional']:
            if TF in uncond_TF_GRE_genes:
                uncond_TF_GRE_genes[TF].update(set(cond_GRE_genes['unconditional'][GRE]))
            else:
                uncond_TF_GRE_genes[TF] = set(cond_GRE_genes['unconditional'][GRE])


all_cond_genes = set.union(*cond_TF_GRE_genes.values())
all_uncond_genes = set.union(*uncond_TF_GRE_genes.values())


#create conditional gene-TF mapping 
all_TFs = set(cond_TF_GRE_genes.keys())
cond_gene_TF = {}
for gene in all_cond_genes:
    for TF in all_TFs:
        if gene in cond_TF_GRE_genes[TF]:
            if gene in cond_gene_TF:
                cond_gene_TF[gene].update({TF})
            else:
                cond_gene_TF[gene] = {TF}

#create gene_TF mapping 
all_TFs = set(uncond_TF_GRE_genes.keys())
gene_TF = {}
for gene in all_uncond_genes:
    for TF in all_TFs:
        if gene in uncond_TF_GRE_genes[TF]:
            if gene in gene_TF:
                gene_TF[gene].update({TF})
            else:
                gene_TF[gene] = {TF}
                


#create modules with hypoxia specificity
cond_modules = {}
for gene, TF_list in cond_gene_TF.items():
    TFs = tuple(TF_list)
    if TFs in cond_modules:
        cond_modules[TFs].append(gene)
    else:
        cond_modules[TFs] = [gene]
        
        
#create modules with MAST GRE-gene mapping
uncond_modules = {}
for gene, TF_list in gene_TF.items():
    TFs = tuple(TF_list)
    if TFs in uncond_modules:
        uncond_modules[TFs].append(gene)
    else:
        uncond_modules[TFs] = [gene]


        
uncond_module_regulator_sizes = [len(x) for x in uncond_modules.keys()]
cond_module_regulator_sizes = [len(x) for x in cond_modules.keys()]
uncond_module_sizes = [len(x) for x in uncond_modules.values()]
cond_module_sizes = [len(x) for x in cond_modules.values()]    




fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
#ax = fig.add_subplot(3,1,1)
#ax.hist(module_regulator_sizes, bins = 57)
#ax.set_xlabel('# TFs/module')
#ax.set_ylabel('Frequency')
#ax.set_title('Unconditional regulators ({})'.format(len(set.union(*[set(x) for x in uncond_modules.keys()]))))
#ax = fig.add_subplot(3,1,2)
#ax.hist(module_sizes, bins = 44)
#ax.set_xlabel('# genes/module')
#ax.set_ylabel('Frequency')
#ax.set_title('Unconditional targets ({})'.format(len(set.union(*[set(x) for x in uncond_modules.values()]))))
ax =  fig.add_subplot(1,1,1)
ax.scatter(cond_module_sizes, cond_module_regulator_sizes)
ax.set_xlabel('# genes/module')
ax.set_ylabel('# TFs/module')
ax.set_title('{} modules'.format(len(cond_module_sizes)))
#ax = fig.add_subplot(3,2,2)
#ax.hist(module_hypox_regulator_sizes, bins = 57)
#ax.set_xlabel('# TFs/module')
##ax.set_ylabel('Frequency')
#ax.set_title('Hypoxia regulators ({})'.format(len(set.union(*[set(x) for x in cond_modules.keys()]))))
#ax = fig.add_subplot(3,2,4)
#ax.hist(module_hypox_sizes, bins = 80)
#ax.set_xlabel('# genes/module')
##ax.set_ylabel('Frequency')
#ax.set_title('Hypoxia targets ({})'.format(len(set.union(*[set(x) for x in cond_modules.values()]))))
#ax =  fig.add_subplot(3,2,6)
#ax.scatter(module_hypox_sizes, module_hypox_regulator_sizes)
#ax.set_xlabel('# genes/module')
#ax.set_ylabel('# TFs/module')
#ax.set_title('{} modules'.format(len(module_hypox_sizes)))
#plt.tight_layout()
plt.savefig('module_sizes_allTFs_matches.eps')
    

#TF-connections per gene

gene_TF_conn_count = {}
for g, t in cond_gene_TF.items():
    gene_TF_conn_count[g] = len(t)
TF_gene_conn_count = {}
for t, g in cond_TF_GRE_genes.items():
    TF_gene_conn_count[t] = len(g)

fig = plt.figure(figsize = (11.69,8.27), dpi = 300)
ax = fig.add_subplot(211)
ax.hist(gene_TF_conn_count.values(), bins = max(gene_TF_conn_count.values()))
ax.set_title("TFs per gene")
ax.set_ylabel('frequency')
ax.vlines(np.mean(list(gene_TF_conn_count.values())),0,ax.get_ylim()[1], color = 'red', label = 'mean')
ax.legend()

#gene-connections per TF           

ax = fig.add_subplot(212)
ax.hist(TF_gene_conn_count.values(), bins = max(TF_gene_conn_count.values()))
ax.vlines(np.mean(list(TF_gene_conn_count.values())),0,ax.get_ylim()[1], color = 'red', label = 'mean')
ax.set_title("genes per TF")
ax.set_ylabel('frequency')
ax.legend()
plt.savefig('TFs_per_gene_and_vice_versa_all_TFs_hypoxia.eps')

#to compare: same spiel for TF binding data

chip_dict = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_dict.p', 'br'))
chip_dict = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/TF_diffex_dic.p', 'br'))

all_default_genes = set.union(set([x for y in chip_dict.values() for x in y]))

all_TFs = set(chip_dict.keys())
gene_TF = {}
for gene in all_default_genes:
    for TF in all_TFs:
        if gene in chip_dict[TF]:
            if gene in gene_TF:
                gene_TF[gene].update({TF})
            else:
                gene_TF[gene] = {TF}

modules = {}
for gene, TF_list in gene_TF.items():
    TFs = tuple(TF_list)
    if TFs in modules:
        modules[TFs].append(gene)
    else:
        modules[TFs] = [gene]

module_regulator_sizes = [len(x) for x in modules.keys()]
module_sizes = [len(x) for x in modules.values()] 

fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
ax =  fig.add_subplot(1,1,1)
ax.scatter(module_sizes, module_regulator_sizes)
ax.set_xlabel('# genes/module')
ax.set_ylabel('# TFs/module')
ax.set_title('{} modules'.format(len(module_sizes)))
plt.savefig('module_sizes_TFOE_bound.eps')
    



gene_TF_conn_count = {}
for g, t in gene_TF.items():
    gene_TF_conn_count[g] = len(t)
TF_gene_conn_count = {}
for t, g in chip_dict.items():
    TF_gene_conn_count[t] = len(g)

fig = plt.figure(figsize = (11.69,8.27), dpi = 300)
ax = fig.add_subplot(211)
ax.hist(gene_TF_conn_count.values(), bins = max(gene_TF_conn_count.values()))
ax.set_title("TFs per gene")
ax.set_ylabel('frequency')
ax.vlines(np.mean(list(gene_TF_conn_count.values())),0,ax.get_ylim()[1], color = 'red', label = 'mean')
ax.legend()         
ax = fig.add_subplot(212)
ax.hist(TF_gene_conn_count.values(), bins = max(TF_gene_conn_count.values()))
ax.vlines(np.mean(list(TF_gene_conn_count.values())),0,ax.get_ylim()[1], color = 'red', label = 'mean')
ax.set_title("genes per TF")
ax.set_ylabel('frequency')
ax.legend()
plt.savefig('TFs_per_gene_and_vice_versa_TFOE_bound.svg')

#corem connections per GRE
conn_counts = [len(x) for x in GRE_cor_dic.values()]
hot_GREs = [k for k,x in GRE_cor_dic.items() if len(x) >45]
fig = plt.figure(figsize = (11.69,8.27), dpi = 300)
ax = fig.add_subplot(111)
ax.hist(conn_counts, bins = 68, label = 'all GREs')
ax.set_xlabel("corems per GRE")
#corem connections per GRE for high confidence GREs
GREs = strenge_hits.GRE_ID.values
conn_counts_hc = [len(GRE_cor_dic[x]) for x in set(GREs)]
ax.hist(conn_counts_hc, bins = 68, color = 'orange', label = 'high confidence GREs')
ax.legend()
plt.savefig('corems_per_GRE.svg')
hot_GREs = GREs[[True if x >= 40 else False for x in conn_counts_hc]]

#determine number of associated genes per GRE in condition
condition = 'unconditional'
GRE_cond_gene_count = {}
for GRE in all_GREs:
    if GRE in cond_GRE_genes[condition]:
        GRE_cond_gene_count[GRE] = len(cond_GRE_genes[condition][GRE])
    else:
        GRE_cond_gene_count[GRE] = 0

conn_counts = GRE_cond_gene_count.values()
fig = plt.figure(figsize = (11.69,8.27), dpi = 300)
ax = fig.add_subplot(111)
ax.hist(conn_counts, bins = max(conn_counts), label = 'all GREs')
ax.set_xlabel("genes per GRE in hypoxia")
#corem connections per GRE for high confidence GREs
GREs = strenge_hits.GRE_ID.values
conn_counts_hc = [x for k,x in GRE_cond_gene_count.items() if k in GREs]
ax.hist(conn_counts_hc, max(conn_counts), color = 'orange', label = 'high confidence GREs')
ax.legend()
plt.savefig('genes_per_GRE_in_hypoxia.svg')


#genes_through_corem_0 = set()
#for corem in GRE_cor_dic[0]:
#    genes_through_corem_0.update(set(corems_genes[corem])) #655 gene sind mit GRE0 assoziiert
#genes_through_corem_15 = set()
#for corem in GRE_cor_dic[15]:
#    genes_through_corem_15.update(set(corems_genes[corem])) #496 gene sind mit GRE0 assoziiert

for _,hit in hits.loc[hits.TF_ID == 'Rv2359'].iterrows():
    TF = hit.TF_ID
    GRE_ID = hit.GRE_ID
    print(GRE_ID)
    peak_set = hit.peaks
    TF_motif_shape = hit.TF_shape
    TF_motif_ID = hit.TF_motif_ID
    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF, peak_set = peak_set, TF_motif_shape = TF_motif_shape, TF_motif_ID = TF_motif_ID)



for TF, match in best_TF_maps.items():
    GRE_ID = match.GRE_ID
    print(GRE_ID)
    peak_set = match.peaks
    TF_motif_shape = match.TF_shape
    TF_motif_ID = match.TF_motif_ID
    tomtom_candidates(GRE_ID = GRE_ID, TF_ID = TF, peak_set = peak_set, TF_motif_shape = TF_motif_shape, TF_motif_ID = TF_motif_ID)