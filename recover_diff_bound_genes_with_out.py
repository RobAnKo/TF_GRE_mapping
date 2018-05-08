#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 11:02:05 2018

@author: rkoch
"""
import subprocess
import pickle
import os
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
nan = np.float('NaN')
import scipy.stats as ss
import matplotlib_venn as venn
import matplotlib.pyplot as plt
import pandas as pd
#the sets of genes we hope to recover using FIMO with motifs derived from peaks outside the promoter regions
diffex_bound = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/TF_diffex_bound.p','rb'))
TFs = diffex_bound.keys()






def recovery_diffex_bound_genes(TF_ID, motif_ID, thresh, sett = 'out'):
      
    
    
    # run analysis against all promoters
    bgfile = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters_bgmodel.txt'
    database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters.fna'  
    
    genes_chip = set(diffex_bound[TF])
    no_diffex_bound = len(genes_chip)
    if no_diffex_bound == 0:
        maxno = 100
    else:
        maxno = 2*no_diffex_bound
    TF_file = '/Users/rkoch/Documents/Data_and_Scripts/MEME_OUTPUT/{}/{}_nonpal/meme.txt'.format(sett, TF_ID)
    if not os.path.isfile(TF_file):
        print('TF {} does not have a motif file'.format(TF_ID))
        return [TF_ID, motif_ID, no_diffex_bound, 0, 0, 1, nan, nan, thresh, 0]
        
    
    
    with open(TF_file, 'r') as TF_f:
        TF_data = TF_f.readlines()
    if TF_data == []:
        print('TF {} has an empty motif file in set {}'.format(TF_ID, sett))
        return [TF_ID, motif_ID, no_diffex_bound, 0, 0, 1, nan, nan, thresh, 1]
    loi = [x for x in TF_data if x.startswith('MOTIF  {}'.format(motif_ID))]
    loi = loi[0].split()
    voi = [float(loi[i]) for i in [4,-1]]
    TF_command =  ['/Users/rkoch/Programs/meme/bin/fimo', 
                   '--bgfile', bgfile, 
                   '--oc', 'temp_out', 
                   '--parse-genomic-coord',
                   '--thresh', str(thresh), 
                   '--max-stored-scores', str(maxno),
                   '--motif', str(motif_ID),
                   '--verbosity', str(1),   
                   TF_file,
                   database]
    sp = subprocess.run(TF_command)
    if sp.returncode != 0:
        print("FIMO run failed for TF %s, motif ID %s"%(TF_ID, motif_ID))
        return [TF_ID, motif_ID, no_diffex_bound, 0, 0, 1, voi[0], voi[1], thresh, 2]
    
    with open('temp_out/fimo.txt') as TF_out:
        TF_hits = TF_out.readlines()[1:]
        TF_hits = [x.split('\t') for x in TF_hits]
    if TF_hits == []:
        print('no significant matches were found for TF %s motif ID %s'%(TF_ID, motif_ID))
        no_TF_hits = 0
        return [TF_ID, motif_ID, no_diffex_bound, 0, 0, 1, voi[0], voi[1], thresh,3]
    else:
        genes_FIMO = set([x[2] for x in TF_hits])
        no_TF_hits = len(genes_FIMO)
        
        no_total = 3996
        no_both = len(genes_FIMO.intersection(genes_chip))
        hpd = ss.hypergeom(no_total, no_TF_hits, no_diffex_bound)
        pv = hpd.sf(no_both-1)
        
           #TF ID, motif_ID, number of genes found to be diffex and bound, no of FIMO gene hits, number of overlapping genes, pvalue of hypergeometric test, length of motif, evalue of motif, threshold used to define 
           #significant FIMO hits, case of analysis(0 or 1: no TF_motif given, 2: TF_motif found, but FIMO failed, 3: TF_motif found, but no  significant hits, 4: TF_motif found, hits found, all good)
    return [TF_ID, motif_ID, no_diffex_bound, no_TF_hits, no_both, pv, voi[0], voi[1], thresh, 4]







#TFs_wo_motif = set(result.loc[result.case_analysis.isin([0,1])].TF_ID)
#TFs_wo_sig_hits = set(result.loc[result.case_analysis.isin([3])].TF_ID)
#TFs_w_significant_recovery = set(result.loc[result.pv <= 0.05].TF_ID)
#venn.venn3(subsets = [TFs_wo_motif,TFs_wo_sig_hits, TFs_w_significant_recovery], set_labels = ['TFs_wo_motif','TFs_wo_sig_hits', 'TFs_w_significant_recovery'])

with PdfPages('pvalues_of_diffex_bound_recovery.pdf') as pp:
    fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
    j = 1
    for sett in ['out', 'all', 'in','in_and_expressed']:
        out = [nan]*len(TFs)*2
        i = 0
        for TF in TFs:
            for motif_id in [1,2]:
                out[i] = recovery_diffex_bound_genes(TF, motif_id, 0.95, sett = sett)
                i += 1
        cols = ['TF_ID', 'motif_ID', 'no_diffex_bound', 'no_TF_hits', 'no_both', 'pv', 'length_TF_motif', 'evalue_TF_motif', 'sig_thresh', 'case_analysis']
        result = pd.DataFrame(out, columns = cols)
        result.to_csv(open('recovery_diff_bound_genes_with_{}_motifs.csv'.format(sett),'w'))
        ax = fig.add_subplot(2,2,j)
        ax.hist(result.pv, bins= np.logspace(np.log10(min(result.pv)), np.log10(max(result.pv)), 100))
        ax.set_xscale('log')
        ax.set_title('Distribution of pvalues when recovering \n differentially expressed and bound genes using FIMO \n on {} motifs'.format(sett))
        ax.set_xlabel('pvalues')
        ax.set_ylabel('frequency')
        j+=1
    plt.tight_layout()
    pp.savefig()
