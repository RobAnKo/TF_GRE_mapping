#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 14:06:26 2018

@author: rkoch

produces tomtom result in html format for a pair of GRE and TF motifs


"""
import subprocess
import os
def tomtom_candidates(GRE_ID = 0, peak_set = 'all', TF_ID = 'Rv0022c', TF_motif_shape = 'pal', TF_motif_ID = 1, e_thresh = 1, dist_method = 'pearson', min_overlap = 5):
    """
    This function takes interesting candidate pairs from matching analysis and produces graphical output for a rerun of 
    those pairs
    """
    
    
    
    cwd = '/Users/rkoch/Documents/Data_and_Scripts/'
    GRE_file = '{:04.0f}_memeOut.txt'.format(GRE_ID) 
    GRE_path = cwd + 'motif_clusters_24/' + GRE_file
    if peak_set not in ['all', 'out', 'in', 'in_and_expressed', 'in_ex']:
        print("Error: peak_set has to be one of 'all', 'out', 'in', 'in_and_expressed' or the short form 'in_ex'")
        return
    if peak_set == "in_ex":
        peak_set = 'in_and_expressed'
    if not TF_ID.startswith('Rv'):
        print("Error: TF_ID has to be a valid TF ID, starting with 'Rv' followed by 4 digits and an optional specifier")
        return
    if not TF_motif_shape in ['pal', 'nonpal']:
        print("Error: TF_motif_shape has to be one of 'pal' or 'nonpal'")
        return
    if not os.path.exists(GRE_path):
        print("Error: the file {} does not exist!".format(GRE_path))
        return
    TF_path = cwd + 'MEME_OUTPUT/' + peak_set + '/%s_%s/meme.txt'%(TF_ID, TF_motif_shape)
    if not os.path.exists(TF_path):
        print("Error: the file {} does not exist!".format(TF_path))
        return    
    
    
    
        
    if TF_motif_ID in [1,2]:
        command = ['/Users/rkoch/Programs/meme/bin/tomtom',
                   '-no-ssc', 
                   '-oc', '.',
                   '-verbosity', '5',
                   '-evalue',
                   '-thresh', '%.3f' % e_thresh,        
                   '-dist', dist_method,
                   '-m', str(TF_motif_ID),
                   '-min-overlap', '%s'%min_overlap,
                   TF_path, GRE_path]
    elif TF_motif_ID == 'both':
        command = ['/Users/rkoch/Programs/meme/bin/tomtom',
                   '-no-ssc', 
                   '-oc', '.',
                   '-verbosity', '5',
                   '-evalue',
                   '-thresh', '%.3f' % e_thresh,        
                   '-dist', dist_method,
                   '-min-overlap', '%s'%min_overlap,
                   TF_path, GRE_path]
    #collect tomtom output
    subprocess.run(command)
    os.rename('tomtom.html', '%s_vs_%s_%s_%s_%s.html'%(GRE_ID, TF_ID, peak_set, TF_motif_shape, TF_motif_ID))
    return

#good_cands = df containing interesting candidates
#cands = lists['all']#.loc[lists['all'].case_of_overlap == 3]
#cands.sort_values(by= 'ratio_rank', inplace=True)
#good_cands = cands.iloc[0:20]
#for values in good_cands.iterrows():
#    vals = values[1].tolist()
#    tomtom_candidates(GRE_ID=vals[0],peak_set='all', TF_ID=vals[1], TF_motif_shape=vals[2], TF_motif_ID= vals[3])