##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Thu Feb 15 13:46:13 2018
#
#@author: rkoch
#"""
#
import os
import pickle
from math import * 
import matplotlib
import matplotlib.pyplot as plt
from statistics import mean 
import progressbar as pb
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sbn
from matplotlib.backends.backend_pdf import PdfPages


plt.rcParams.update({'font.size': 12})


"""
#print(os.getcwd())
#
#
#
#
####---------------------------------------------------------------------------------------------###
#                    #check how many genes are used to create GRE using MEME
#                    # and why are they created using genes anyway? Aren't they just averages of CRMs?
#                    # also: what is the length of these genes (might be promoter region?)
#                    
os.chdir('/Users/rkoch/Documents/Data_and_Scripts/motif_clusters_24')

l = os.listdir()
l = [x for x in l if x.endswith('memeOut.txt')]
nos_of_genes = []
len_of_genes = []
minws = []
maxws = []
n = 0
for file in l:
    with  open(file) as f_acc:
        lines = f_acc.readlines()
        line1 = [x for x in lines if x.startswith('data:')] #line containing number of genes
        no = line1[0].split(" ")[-1].split('\n')[0]
        nos_of_genes.append(int(no))
        line2 = [x for x in lines if x.startswith("Sequence name            Weight")][0] #beginning of list of genes
        start_ind = lines.index(line2) + 2
        end_ind = start_ind + ceil(float(no)/2)
        for ind in range(start_ind,end_ind):
            line3 = lines[ind]
            line3 = line3.split(' ')
            line3 = [x for x in line3 if x != '']
            inds = [2,5] # where the numbers live
            if len(line3) < 6:
                lens_of_genes = [int(line3[2])]
            else:
                lens_of_genes = [int(line3[x]) for x in inds]
            len_of_genes.extend(lens_of_genes)
        n= n+1
        line4 = [x for x in lines if x.startswith("command")][0].split('-minw ')[1].split(' ')
        minw, maxw = int(line4[0]), int(line4[2])
        minws.append(minw)
        maxws.append(maxw)
        
#print(min(minws), max(minws))
#print(min(maxws), max(maxws))
plt.hist(nos_of_genes, bins=50)
plt.title('Numbers of genes used to create GRE PSSMs\n Total #GREs: %s\n Total #Genes %s' %(n, sum(nos_of_genes)))
plt.ylabel('frequency')
plt.xlabel('number of genes')
plt.show()


plt.hist(len_of_genes, bins = 35)
plt.title('Length of all genes used to create GRE PSSMs\n Total #genes: %s'%sum(nos_of_genes))
plt.ylabel('frequency')
plt.xlabel('length of gene')
plt.show()
#
####---------------------------------------------------------------------------------------------###
#                    #check how long the GREs are
#                    
GRE_lengths = {}                    
for file in l:
    GRE_ID = int(file.split('_')[0])
    with  open(file) as f_acc:
        lines = f_acc.readlines()
        line = [x for x in lines if x.startswith('MOTIF  1 MEME')][0] #line containing width of motif
        line = line.split(' ')
        empty = line.count('')
        for i in range(empty):
            line.remove('')
        no = int(line[4])
        GRE_lengths[GRE_ID] = no
#        
#os.chdir('/Users/rkoch/Documents/Data_and_Scripts/')
#
####---------------------------------------------------------------------------------------------###
#                    #check the TF PSSM stuff
#
#os.chdir("/Users/rkoch/Documents/Data_and_Scripts/MEME_output/out")
#
#l = os.listdir() #a list of all TF PSSM files
#l = [x for x in l if x.startswith('Rv')]
#TFs = [x.split('_')[0] for x in l]
#TFs = set(TFs)
#flavours = ['pal','nonpal']
#nos_of_genes = []
#empty_files = 0
#for TF in TFs:
#    for flavour in flavours:
#        file = '%s_%s/meme.txt'%(TF,flavour)
#        with open(file) as f:
#            lines = f.readlines()
#            if len(lines) > 0:
#                line1 = [x for x in lines if x.startswith('data:')] #line containing number of genes
#                no = line1[0].split(" ")[-1].split('\n')[0]
#                nos_of_genes.append(int(no))
#            else:
#                empty_files = empty_files + 1
#
#plt.hist([x for x in nos_of_genes], bins = 100)
#plt.title('Numbers of genes used to create TF PSSMs\n Total #TFs: %s'%len(TFs))
#plt.ylabel('frequency')
#plt.xlabel('number of genes')
#plt.show()
# 
#
#
####---------------------------------------------------------------------------------------------###
#                    #check how long the TF motifs are
#                    
TF_lengths = {}
os.chdir("/Users/rkoch/Documents/Data_and_Scripts/MEME_output")
sets = ['all', 'out', 'in', 'in_and_expressed']
for set_ in sets:
    TF_lengths[set_] = {'pal':{},'nonpal':{}}
    l = [x for x in os.listdir(set_) if 'Rv' in x] #a list of all TF PSSM files
    for TF_ in l:
        spl = TF_.split('_')
        TF = spl[0]
        p = spl[1]
        with open('%s/%s/meme.txt'%(set_,TF_)) as f_acc:
            lines = f_acc.readlines()
            lines = [x for x in lines if x.startswith('letter-probability matrix')] #line containing width of motif
            for line in lines:
                line = line.split(' ')
                empty = line.count('')
                for i in range(empty):
                    line.remove('')
                no = int(line[5])
                TF_lengths[set_][p][TF] = no

TF_pal_lengths = [x for set_ in TF_lengths for style in TF_lengths[set_] for x in TF_lengths[set_][style].values() if style =='pal']
TF_nonpal_lengths = [x for set_ in TF_lengths for style in TF_lengths[set_] for x in TF_lengths[set_][style].values() if style =='nonpal']
plt.hist([list(GRE_lengths.values()), TF_pal_lengths, TF_nonpal_lengths], bins = 30, histtype="bar")
plt.legend(['GRE', 'TF_pal', 'TF_nonpal'])
plt.title('lengths of motifs')
plt.xlabel('length [nt]')
plt.ylabel('frequency')
plt.tight_layout()






        
os.chdir('/Users/rkoch/Documents/Data_and_Scripts')                   
                    

###---------------------------------------------------------------------------------------------###

#check the lengths and number and evals of  motifs
os.chdir('/Users/rkoch/Documents/Data_and_Scripts/MEME_OUTPUT')                   
#os.chdir('/home/karlchen/ISB_data/MEME_OUTPUT')
types = ['all', 'in', 'in_and_expressed', 'out']
palindromic = ['_pal', '_nonpal']
lengths = {}               # to compare lengths of motifs between different sets (pal,nonpal,all,in,etc)
numbers = {}               # to count how many TFs per set (in, out etc)
evals = {'1_pal':[], '2_pal': [], '1_nonpal':[], '2_nonpal': []}  # to compare the evalues of motif IDs
#evals_pal_vs_nonpal = {'_pal':[], '_nonpal'[]}
for typ in types:
    for p in palindromic:
        lengths['%s%s'%(typ, p)] = []
        numbers['%s%s'%(typ, p)] = 0
        TF_files = [typ + '/' + x + '/meme.txt' for x in os.listdir(typ) if x.endswith(p)] 
        for TF_file in TF_files:
            with open(TF_file, 'r') as f:
                content = f.readlines()
                if not content:
                    print('this file is empty: %s'%TF_file)
                else:
                    data_lines = [line for line in content if line.startswith('MOTIF')]
                    chunks =  [line.split() for line in data_lines]
                    important_chunks = [(x[1], float(x[-1]), int(x[4])) for x in chunks] #(motifID, evalue, length)
                    #print(len(data_lines))
                    numbers['%s%s'%(typ, p)] += len(data_lines)
                    for x in important_chunks:
                        evals['%s%s'%(x[0], p)].append(x[1] + np.nextafter(0,1))
                        lengths['%s%s'%(typ, p)].append(x[-1])
#                        
#
#
#
labels = []
data = []
for k, v in lengths.items():
    labels.append(k)    #different sets
    data.append(v)      # corresponding length distributions   
labels = [x for _, x in sorted(zip(data,labels), key=lambda pair: mean(pair[0]))] # sort sets based on mean length of motifs in set
new_labels = []
for l in labels:
    new_labels.append(l + ': '+ str(numbers[l]))
data = sorted(data,key = lambda x: mean(x))
new_labels = [0] + ['GREs: 152'] + new_labels
data = [GRE_lengths] + data



with PdfPages('sdfTF_motif_lengths.pdf') as pdfz:
    plt.figure(figsize=(11.69, 8.27), dpi=100)
    ax = plt.subplot(111)
    ax.set_yticks(range(0,10))
    ax.set_yticklabels(new_labels)
    plt.violinplot(data, vert=False, showextrema=False, showmeans=True)
    plt.xlabel('lengths of TF motifs [nt]')
    plt.ylabel('different motif sets')
    plt.title('Distribution of TF motif length for different motif sets')
    plt.text(s='set:number of motifs', x = .1, y = 8.5)
    plt.tight_layout()
    pdfz.savefig()
#
#
#
# plot log10 of evalues
log_evals = {'1_pal':[], '2_pal':[], '1_nonpal':[], '2_nonpal':[]}
for k, v in evals.items():
    log_evals[k] = np.log10(v)
    
with PdfPages('Evalues_pal_vs_nonpal_vs_1_vs_2.pdf') as pdfy:
    fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
    ax = fig.add_subplot(111)
    plt.violinplot(list(log_evals.values()), showmeans= True, showextrema=False, showmedians=False)
    plt.title('Distribution of E-values for different motif sets')
    plt.xlabel('Set')
    plt.ylabel('log10(E-value)')
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(list(log_evals.keys()))
    plt.tight_layout()
    pdfy.savefig()
    
    fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_ylim([-3,10])
    plt.violinplot(list(log_evals.values()), showmeans= True, showextrema=False, showmedians=False)
    plt.title('Distribution of E-values for different motif sets')
    plt.xlabel('Set')
    plt.ylabel('log10(E-value)')
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(list(log_evals.keys()))
    plt.tight_layout()
    pdfy.savefig()
"""    
 ###---------------------------------------------------------------------------------------------###
"""
#analyse the matchlists (motif-based)
#os.chdir('/home/karlchen/ISB_data/Dictionaries')
os.chdir('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries')
ml_all = pickle.load(open('combined_match_list_with_center_evals_noseqs_all.p', 'rb'))
ml_in = pickle.load(open('combined_match_list_with_center_evals_noseqs_in.p', 'rb'))
ml_out = pickle.load(open('combined_match_list_with_center_evals_noseqs_out.p', 'rb'))
ml_in_and_expressed = pickle.load(open('combined_match_list_with_center_evals_noseqs_in_and_expressed.p', 'rb'))

#ml_all = pickle.load(open('combined_match_list_all.p', 'rb'))
#ml_in = pickle.load(open('combined_match_list_in.p', 'rb'))
#ml_out = pickle.load(open('combined_match_list_out.p', 'rb'))
#ml_in_and_expressed = pickle.load(open('combined_match_list_in_and_expressed.p', 'rb'))

lists = {'in':ml_in, 'in_ex':ml_in_and_expressed, 'out':ml_out, 'all': ml_all}
cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'case_of_overlap', 'length_overlap', 'length_GRE_uncut', 'length_GRE_cut', 'number_spliced_nts_GRE', 'length_TF_uncut', 'length_TF_cut', 'number_spliced_nts_TF', 'pn', 'pv', 'orientation', 'center', 'evalue_TF', 'no_seqs_TF', 'evalue_GRE', 'no_seqs_GRE']
short_cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID','case_of_overlap', 'length_overlap','number_spliced_nts_GRE', 'length_GRE_uncut', 'length_GRE_cut', 'length_TF_uncut', 'pn', 'pv', 'orientation']
#store all the results as dataframe
for k in lists.keys():
    lists[k] = pd.DataFrame(lists[k], columns = cols)
    lists[k].to_csv(open('motif_matches_{}.csv'.format(k), 'w'))


for key in lists.keys():
    #optional: reduce matches to those with motifID 1
    #lists[key] = lists[key].loc[lists[key].TF_motif_ID == 1] 
    #lists[key].index = range(len(lists[key]))
    
    #'improvement' of pvalue
    lists[key]['pv_ratio_GRE_cut'] =  lists[key].pn/lists[key].pv   
    
    #rank based on normal pvalue
    lists[key].sort_values(by = "pn", inplace = True)
    lists[key]['pn_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pn_rank', index = lists[key].index)  
    #rank based on cut pvalue
    lists[key].sort_values(by = 'pv', inplace = True)
    lists[key]['pv_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pv_rank', index = lists[key].index)  
    #rank based on pv_ratio_GRE_cut
    lists[key].sort_values(by = 'pv_ratio_GRE_cut', ascending = False, inplace = True)
    lists[key]['pv_ratio_GRE_cut_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pv_ratio_GRE_cut_rank', index = lists[key].index)
    
    #lists[key].sort_values(by= 'pv_rank', inplace = True)
    #lists[key].set_index(['GRE_ID', 'TF_ID', 'TF_shape'], inplace=True)
    lists[key].sort_index(inplace = True)



sublists_GRE = {} # a dictionary: for all GREs a list of the matched TFs
sublists_TF = {} # a dictionary: for all TFs a list of the matched GREs
TF_IDs = set(lists['all'].TF_ID)
GRE_IDs = set(lists['all'].GRE_ID)
progbar = pb.ProgressBar()
for key in progbar(lists.keys()):
    print(key)
    sublists_GRE[key] = {} 
    sublists_TF[key] = {}
    for GI in GRE_IDs:
        sublists_GRE[key][GI] = lists[key].loc[lists[key].GRE_ID == GI]
        sublists_GRE[key][GI].sort_values(by = 'pv', inplace = True)
    for TF in TF_IDs:
        print(TF)
        sublists_TF[key][TF] = lists[key].loc[lists[key].TF_ID == TF]
        sublists_TF[key][TF].sort_values(by = 'pv', inplace = True)


###should i subdivide for pal\nonpal\1\2\?
        
#reciprocal matching
        
reciproc_matches = {}
for key in sublists_TF.keys():
    reciproc_matches[key] = {}
    progling = pb.ProgressBar()
    for n in progling(range(1,5)):
        reciproc_matches[key][n] = []
        top_TFs = {}
        top_GREs = {}
        for GRE in GRE_IDs:
            top_TFs[GRE] = sublists_GRE[key][GRE].iloc[0:n]
        for TF in TF_IDs:
            top_GREs[TF] = sublists_TF[key][TF].iloc[0:n]
        for G, Ts in top_TFs.items():
            for T in Ts.TF_ID:
                indices = np.where(top_GREs[T].GRE_ID.values == G)[0]
                for i in indices:
                    TF_shape = Ts.TF_shape.iloc[i]
                    TF_motif_ID = Ts.TF_motif_ID.iloc[i]
                    reciproc_matches[key][n].append((G, key, T, TF_shape, TF_motif_ID))
                    print('yaay')

reci_GREs = {}
all_GREs = set(lists['all'].GRE_ID)
reci_TFs = {}
all_TFs = set(lists['all'].TF_ID)
ratioGRE = {}
ratioTFs = {}
for key in reciproc_matches.keys():
    reci_GREs[key] = {}
    reci_TFs[key] = {}
    ratioGRE[key] = {}
    ratioTFs[key] = {}
    for level in range(1,5):
        matches = reciproc_matches[key][level]    
        reci_GREs[key][level] = set([x[0] for x in matches])
        reci_TFs[key][level] = set([x[2] for x in matches])
        ratioGRE[key][level] = len(reci_GREs[key][level])/len(all_GREs)
        ratioTFs[key][level] = len(reci_TFs[key][level])/len(all_TFs)


    
#unique_reciprocal_matches = set([x for k in reciproc_matches.keys() for n in reciproc_matches[k].keys() for x in reciproc_matches[k][n]])
#reciproc_GREs = [x[0] for x in unique_reciprocal_matches]
#reciproc_TFs = [x[1] for x in unique_reciprocal_matches]

#len_reciproc_GREs =  
#len_reciproc_TFs = 








pp =  PdfPages('pvalue_ratio_vs_cut_nucleotides.pdf')
i = 1
fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    df_filt = df.loc[df.pn<0.01]
    ax = fig.add_subplot(2,2,i)
    ax.plot(df.number_spliced_nts_GRE, df.pv_ratio_GRE_cut, 'x', markersize = 0.5, label = 'all matches', color = '#deebf7')
    ax.plot(df_filt.number_spliced_nts_GRE, df_filt.pv_ratio_GRE_cut, 'x', markersize = 1, label = 'initial pvalue <0.01', color = '#3182bd')
    ax.hlines(1, 0, 25, label = 'no change')
    ax.set_yscale('log')
    ax.set_ylabel('pv_ratio_GRE_cut: original pvalue/pvalue after splicing')
    ax.set_xlabel('# of nucleotides cut off')
    plt.legend()
    ax.set_title('Change in tomtom p-values for TF-GRE matches \n motif set: %s'%key)
    plt.tight_layout()
    i += 1
pp.savefig()
pp.close()

#plot pv_ratio_GRE_cuts vs overlap

pp =  PdfPages('overlap_vs_pvalues.pdf')
i = 1
fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    df_filt = df.loc[df.pn<0.01]
    ax = fig.add_subplot(2,2,i)
    ax.plot(df.length_overlap, df.pn, 'x', markersize = 0.5, label = 'all matches', color = '#deebf7')
    ax.plot(df_filt.length_overlap, df_filt.pn, 'xr', markersize = 1, label = 'initial pvalue <0.01', color = '#3182bd')
    ax.set_yscale('log')
    plt.ylabel('original pvalue')
    plt.xlabel('# of nucleotides overlap')
    plt.legend()
    plt.title('tomtom p-values for TF-GRE matches with different overlaps \n motif set: %s'%key)
    plt.tight_layout()
    i += 1
pp.savefig()
pp.close()

#plot pv_ratio_GRE_cuts vs coverage

pp =  PdfPages('coverage_vs_pv_ratio_GRE_cuts.pdf')
i = 1
fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    df_filt = df.loc[df.pn<0.01]
    ax = fig.add_subplot(2,2,i)
    ax.plot(df.TF_coverage, df.pv_ratio_GRE_cut, 'x', markersize = 0.5, label = 'all matches', color = '#deebf7')
    ax.plot(df_filt.TF_coverage, df_filt.pv_ratio_GRE_cut, 'x', markersize = 1, label = 'initial pvalue <0.01', color = '#3182bd')
    ax.set_yscale('log')
    plt.ylabel('pv-ratio: $p_{normal}$/$p_{cut}$')
    plt.xlabel('fraction of TF bound')
    plt.legend()
    plt.title('tomtom p-value changes for TF-GRE matches \n vs TF coverage \n motif set: %s'%key)
    plt.subplots_adjust(left  = 0.125, right = 0.9, bottom = 0.1, top = 0.9, wspace = 0.2, hspace = 0.2)
    #plt.tight_layout()
    i += 1
pp.savefig()
pp.close()


#compare palindromic and nonpalindromic matches
pp =  PdfPages('pvalues_pal_vs_nonpal.pdf')
f = plt.figure(figsize=(11.69, 8.27), dpi=300)
i = 1
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    ax = f.add_subplot(2,2,i)
    v1 = df.pn.loc[df.TF_shape == 'pal']
    v2 = df.pn.loc[df.TF_shape == 'nonpal']  
    ax.violinplot([v1, v2],vert=False, showextrema=False, showmeans=True)
    ax.set_yticks([1,2])
    ax.set_yticklabels(['pal', 'nonpal'])
    plt.xlabel('pvalue')
    ax.set_xscale('log')
    ax.set_title('original pvalues for palindromic and nonpalindromic motif matches \n for set: %s'%key)
    plt.plot()
    i+=1
    plt.tight_layout()
pp.savefig()
pp.close()

pp =  PdfPages('new_pvalues_pal_vs_nonpal.pdf')
f = plt.figure(figsize=(11.69, 8.27), dpi=300)
i = 1
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    ax = f.add_subplot(2,2,i)
    v1 = df.pv.loc[df.TF_shape == 'pal']
    v2 = df.pv.loc[df.TF_shape == 'nonpal']  
    ax.violinplot([v1, v2],vert=False, showextrema=False, showmeans=True)
    ax.set_yticks([1,2])
    ax.set_yticklabels(['pal', 'nonpal'])
    plt.xlabel('pvalue')
    ax.set_xscale('log')
    ax.set_title('new pvalues for palindromic and nonpalindromic motif matches \n for set: %s'%key)
    plt.plot()
    i+=1
    plt.tight_layout()
pp.savefig()
pp.close()


pp =  PdfPages('pvalues_different_cases.pdf')
f = plt.figure(figsize=(11.69, 8.27), dpi=300)
i = 1
for key in lists.keys():
    df = lists[key].loc[lists[key].pv_ratio_GRE_cut.notnull()]
    ax = f.add_subplot(2,2,i)
    v1 = df.pn.loc[df.case_of_overlap == 1]
    v2 = df.pn.loc[df.case_of_overlap == 2]
    v3 = df.pn.loc[df.case_of_overlap == 3]
    v4 = df.pn.loc[df.case_of_overlap == 4]
    ax.violinplot([v1, v2, v3, v4],vert=False, showextrema=False, showmeans=True)
    ax.set_yticks([1,2,3,4)
    ax.set_yticklabels(['case 1', 'case 2', 'case 3', 'case 4'])
    plt.xlabel('pvalue')
    ax.set_xscale('log')
    ax.set_title('pvalue distributions for all 4 cases of matching \n for set: %s'%key)
    plt.plot()
    i+=1
    plt.tight_layout()
pp.savefig()
pp.close()


for i in [10,20,50,100,200,500,1000, 2000, 10000]:
        
    pn_set = set(df.index[df.pn_rank<=i])
    pv_set = set(df.index[df.pv_rank<=i])
    pv_ratio_GRE_cut_set = set(df.index[df.pv_ratio_GRE_cut_rank<=i])
    venn.venn3([pn_set, pv_set, pv_ratio_GRE_cut_set], set_labels = ['pn', 'pv', 'rank_pv_ratio_GRE_cut'])
plt.show()


 ###---------------------------------------------------------------------------------------------###

# a little excursion: try to predict pvalue from ohter values

df = lists['all']

feat_cols = ['case_of_overlap','length_overlap','length_GRE_uncut','length_GRE_cut','length_TF_uncut','length_TF_cut', 'number_spliced_nts_GRE','number_spliced_nts_TF','evalue_TF', 'no_seqs_TF', 'evalue_GRE', 'no_seqs_GRE']

X = df.loc[:,feat_cols]
y = df.loc[:,'pv']
from sklearn import linear_model
from sklearn import svm
from sklearn.model_selection import cross_val_score

classifiers = [
    #svm.SVR(),
    linear_model.SGDRegressor(),
    linear_model.BayesianRidge(),
    linear_model.LassoLars(),
    linear_model.ARDRegression(),
    linear_model.PassiveAggressiveRegressor(),
    linear_model.TheilSenRegressor(),
    linear_model.LinearRegression()]
scores_out = {}
for item in classifiers:
    print(item)
    print(item.score)
    clf = item
    scores = cross_val_score(clf, X, y, cv = 10)
    scores_out[str(type(item)).split('.')[-1]] = scores



clf = svm.SVC(kernel='linear', C=1)
scores = cross_val_score(clf, iris.data, iris.target, cv=5)
scores      
###---------------------------------------------------------------------------------------------###







#check rank correlations between sets for every GRE based on normal pvalue, cut pvalue or pvalue ratio
GRE_IDs = set(lists['all'].GRE_ID)
kendall_taus_dict = {}
#check rank correlations between sets for every GRE based on normal pvalue
prog = pb.ProgressBar()
for GI in prog(GRE_IDs):
    kendall_taus_dict[GI] = pd.DataFrame(index = [key for key in sublists_GRE], columns = [key for key in sublists_GRE])
    progger = pb.ProgressBar()
    for k1 in progger(sublists_GRE.keys()):
        for k2 in sublists_GRE.keys():
            print(GI, k1, k2)
            kendall_taus_dict[GI][k1].loc[k2],_ = stats.kendalltau(list(sublists_GRE[k1][GI].pv_ratio_GRE_cut), list(sublists[k2][GI].pv_ratio_GRE_cut), nan_policy = 'omit')
pp = PdfPages('kendall_tau_for_all_GREs_pv_ratio_GRE_cut.pdf')
progolont = pb.ProgressBar()
for GI in progolont(kendall_taus_dict.keys()):
    plt.figure(GI)
    plt.title("Kendall's tau correlation for match ranking based on pvalue ratios \n for different motif sets\n for GRE %s"%GI)
    sbn.heatmap(kendall_taus_dict[GI].astype('float'), annot = True)
    plt.tight_layout()
    pp.savefig()
pp.close()









   

##create match lists based on pvalue threshold 
#sign_lists = {}
#for k in lists.keys():
#    sig = lists[k].pv < 0.05
#    sign_lists[k] = lists[k][sig]
#
##plot pn ranks vs pv ranks in pvalue-filtered matches
#pp = PdfPages('filtered_ranks_pn_vs_pv.pdf')   
#for k in sign_lists.keys():
#    f = plt.figure()
#    ax = sign_lists[k].plot.scatter('pn_rank', 'pv_rank',s= 0.001, title = 'ranks of normal vs cut pvalues in dataset %s'%k)
#    pp.savefig(ax.figure)
#
#pp.close()
#
#pp = PdfPages('filtered_pn_vs_pv.pdf')   
#for k in sign_lists.keys():
#    f = plt.figure()
#    ax = sign_lists[k].plot.scatter('pn', 'pv',s= 0.001, title = 'normal vs cut pvalues in dataset %s'%k, logx = True, logy = True)
#    pp.savefig(ax.figure)
#
#pp.close()
#
#pp = PdfPages('filtered_ratio_rank_vs_pn.pdf')   
#for k in sign_lists.keys():
#    f = plt.figure()
#    ax = sign_lists[k].plot.scatter('pv_ratio_GRE_cut_rank', 'pn',s= 0.001, title = 'pv_ratio_GRE_cut rank vs normal pvalues in dataset %s'%k, logx = False, logy = True)
#    pp.savefig(ax.figure)






###---------------------------------------------------------------------------------------------###            
"""
#analyse the output of position/gene based matching
for level in ['genes', 'positions']:
    match_list = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/match_list_{}.p'.format(level), 'br'))
    cols = ['GRE_ID', 'TF_ID', 'n_GRE', 'n_TF', 'n_both', 'pv']
    df = pd.DataFrame(match_list, columns = cols)
    df.to_csv(open('{}_based_matches.csv'.format(level), 'w')) #to store it
    GRE_IDs = set(df.GRE_ID)
    TF_IDs = set(df.TF_ID)
    pvals = df.pv
    alpha = 0.05
    bf = alpha/len(pvals)   
    exp_hits = df.n_TF
    GRE_FIMO_hits = df.n_GRE
    double_hits = df.n_both
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot(exp_hits, GRE_FIMO_hits, double_hits, '.')
    #ax.set_xlabel('# experimentally determined bound genes')
    #ax.set_ylabel('# FIMO found genes')
    #ax.set_zlabel('# experimentally determined genes found by FIMO')
    with PdfPages('pvalues_{}_based_match.pdf'.format(level)) as pp:
        fig = plt.figure(figsize = (11.69, 8.27), dpi=300)
        ax = fig.add_subplot(111)
        ax.hist(pvals, bins = np.logspace(log10(min(pvals[pvals != min(pvals)])), log10(max(pvals)), 100))
        ax.vlines(alpha,0,20000,'r', linestyles = 'dashed',label = '$p = \\alpha = %s$'%alpha, linewidth = 0.5)
        ax.vlines(bf,0,20000,'r', linestyles = 'dotted',label = '$p = \\alpha_{BF}$ = %3.2e'%bf, linewidth = 0.5)
        ax.set_title('Distribution of hypergeom. pvalues for {} {}-based matches'.format(len(pvals), level))
        ax.set_xlabel('pvalue')
        ax.set_ylabel('Frequency')
        ax.set_xscale('log')
        ax.legend()
        ax.text(x = alpha/10000, y = 5000, s = 'matches(p<$\\alpha$):\n {}'.format(len([x for x in pvals if x<alpha])))
        ax.text(x = bf/10000, y = 2600, s = 'matches(p<$\\alpha_{BF}$):\n %s'%len([x for x in pvals if x<bf]))
        plt.tight_layout()
        pp.savefig()
    
    with PdfPages('pvals_vs_exp_hits_{}.pdf'.format(level)) as pp:
        fig = plt.figure(figsize = (11.69, 8.27), dpi=300)
        ax = fig.add_subplot(111)
        ax.scatter(x = exp_hits, y = pvals, marker = 'x', s = 0.1, c = 'k')
        ax.hlines(alpha, 0, max(exp_hits), linestyles = 'dashed', color = 'r', label = '$p = \\alpha = %s$'%alpha, linewidth = 0.5)
        ax.hlines(bf, 0, max(exp_hits), linestyles = 'dotted', color = 'r', label = '$p = \\alpha_{BF}$ = %3.2e'%bf, linewidth = 0.5)
        ax.set_title('pvalues vs. #experimental hits for {} {}-based matches'.format(len(pvals), level))
        ax.set_xlabel('{} bound in ChipSeq experiment'.format(level))
        ax.set_ylabel('pvalue')
        ax.set_ylim([min(pvals[pvals != min(pvals)]), max(pvals)+0.5])
        ax.set_yscale('log')
        ax.legend()
        #ax.text(x = 0.1, y = 2600, s = 'matches(p<{}): {}'.format(alpha, len([x for x in pvals if x<alpha])))
        #ax.text(x = 0.1, y = 2000, s = 'matches(p<{:3f}): {}'.format(bf, len([x for x in pvals if x<bf])))
        plt.tight_layout()
        pp.savefig()

#significant sets
sig_g = [(x.GRE_ID, x.TF_ID) for _, x in ml_gp.iterrows() if x.pv_genes <= bf]
sig_p = [(x.GRE_ID, x.TF_ID) for _, x in ml_gp.iterrows() if x.pv_positions <= bf]
sig_g = set(sig_g)
sig_p = set(sig_p)  
with PdfPages('venn_sig_hits_position.pdf') as pp:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn2([sig_g, sig_p], set_labels = ['genes','positions'])
    ax.set_title('significant matches based on genes/positions')
    pp.savefig()
"""
###---------------------------------------------------------------------------------------------###            
"""
#compare mapping based on motif with mapping based on genes (with mapping based on positions)
    
match_list_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/match_list_genes.p', 'br'))
match_list_positions = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/match_list_positions.p', 'br'))
cols = ['GRE_ID', 'TF_ID', 'n_GRE', 'n_TF', 'n_both', 'pv']
ml_g = pd.DataFrame(match_list_genes, columns = cols)
ml_p = pd.DataFrame(match_list_positions, columns = cols)
os.chdir('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/')
#for both sequences cut
ml_all = pickle.load(open('combined_match_list_with_center_evals_noseqs_all.p', 'rb'))
ml_in = pickle.load(open('combined_match_list_with_center_evals_noseqs_in.p', 'rb'))
ml_out = pickle.load(open('combined_match_list_with_center_evals_noseqs_out.p', 'rb'))
ml_in_and_expressed = pickle.load(open('combined_match_list_with_center_evals_noseqs_in_and_expressed.p', 'rb'))
#for only GRE cut
ml_all_s = pickle.load(open('combined_match_list_with_center_all.p', 'rb'))
ml_in_s = pickle.load(open('combined_match_list_with_center_in.p', 'rb'))
ml_out_s = pickle.load(open('combined_match_list_with_center_out.p', 'rb'))
ml_in_and_expressed_s = pickle.load(open('combined_match_list_with_center_in_and_expressed.p', 'rb'))

lists = {'in':ml_in, 'in_ex':ml_in_and_expressed, 'out':ml_out, 'all': ml_all}
short_lists = {'in':ml_in_s, 'in_ex':ml_in_and_expressed_s, 'out':ml_out_s, 'all': ml_all_s}
cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'case_of_overlap', 'length_overlap', 'length_GRE_uncut', 'length_GRE_cut', 'number_spliced_nts_GRE', 'length_TF_uncut', 'length_TF_cut', 'number_spliced_nts_TF', 'pn', 'pv', 'orientation', 'center', 'evalue_TF', 'no_seqs_TF', 'evalue_GRE', 'no_seqs_GRE']
short_cols = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'case_of_overlap', 'length_overlap', 'number_spliced_nts_GRE', 'length_GRE_uncut', 'length_GRE_cut', 'length_TF_uncut', 'pn', 'pv', 'orientation', 'center']
for k in lists.keys():
    #lists[k] = pd.DataFrame(lists[k], columns = cols)
    lists[k] = pd.DataFrame(lists[k], columns = cols)
    short_lists[k] = pd.DataFrame(short_lists[k], columns = short_cols)
    #store all the results as dataframe
    lists[k].to_csv(open('motif_matches_{}_big.csv'.format(k), 'w'))
    short_lists[k].to_csv(open('motif_matches_{}_only_GRE.csv'.format(k), 'w'))



for key in lists.keys():
    #optional: reduce matches to those with motifID 1
    #lists[key] = lists[key].loc[lists[key].TF_motif_ID == 1] 
    #lists[key].index = range(len(lists[key]))
    
    #add the pvalues for single cut from short lists (i.e. tomtom runs without cutting TF)
    lists[key] = pd.merge(lists[key], short_lists[key][['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID', 'pv']], on = ['GRE_ID', 'TF_ID', 'TF_shape', 'TF_motif_ID'], suffixes = ['_both_cut', '_GRE_cut'])
    
    #'improvement' of pvalue with cut GRE
    lists[key]['pv_ratio_GRE_cut'] =  lists[key].pn/lists[key].pv_GRE_cut
    #fraction of TF motif that is aligned with GRE motif
    lists[key]['TF_coverage'] = lists[key].length_overlap/lists[key].length_TF_uncut
    #rank based on normal pvalue
    lists[key].sort_values(by = "pn", inplace = True)
    #lists[key]['pn_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pn_rank', index = lists[key].index)  
    #rank based on cut pvalue
    lists[key].sort_values(by = 'pv_GRE_cut', inplace = True)
    #lists[key]['pv_GRE_cut_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pv_GRE_cut_rank', index = lists[key].index)  
    #rank based on pv_ratio_GRE_cut
    lists[key].sort_values(by = 'pv_ratio_GRE_cut', ascending = False, inplace = True)
    #lists[key]['pv_ratio_GRE_cut_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pv_ratio_GRE_cut_rank', index = lists[key].index)
    lists[key].sort_index(inplace = True)
    #'improvement' of pvalue with cut GRE and TF
    lists[key]['pv_ratio_both_cut'] = lists[key].pn/lists[key].pv_both_cut
    lists[key].sort_values(by = 'pv_ratio_both_cut', ascending = False, inplace = True)
    #lists[key]['pv_ratio_both_cut_rank'] = pd.Series(range(1, len(lists[key])+1), name = 'pv_ratio_GRE_cut_rank', index = lists[key].index)
#lists['genes'] = ml_g
#lists['positions'] = ml_p

#plot distributions of evalues, pvalues, no_seqs, length_TF

with PdfPages('Evalues.pdf') as ppe, PdfPages('lengths.pdf') as ppl:
    fige = plt.figure(1, figsize=(23.38, 16.54), dpi=300)
    figl = plt.figure(2, figsize=(23.38, 16.54), dpi=300)
    fig
    i = 0
    for key in lists.keys():
        i += 1
        axe = fige.add_subplot(2,2,i)
        axl = figl.add_subplot(2,2,i)
        axe.hist(lists[key].evalue_TF.loc[~np.isnan(lists[key].evalue_TF)])
        axl.hist(lists[key].length_TF_uncut.loc[~np.isnan(lists[key].length_TF_uncut)])
    ppe.savefig(figure = fige)
    ppl.savefig(figure = figl)





    

#compare ratio of pvalues cutting GREs vs cutting both motifs
    fig = plt.figure(figsize=(23.38, 16.54), dpi=300)
    with PdfPages('pval_ratios_GREcut_vs_bothcut.pdf') as pp:
        i = 1
        for key in ['all','out','in','in_ex']:
            ax = fig.add_subplot(2,2,i)
            df = lists[key]
            df_filt = df.loc[df.pn <= 0.01]
            ax.scatter(df.pv_ratio_GRE_cut, df.pv_ratio_both_cut, s = 0.2,color = '#deebf7')
            ax.scatter(df_filt.pv_ratio_GRE_cut, df_filt.pv_ratio_both_cut, s = 0.2,color = '#3182bd')
            ax.set_xlabel('GRE cut')
            ax.set_ylabel('GRE & TF cut')
            ax.set_title('Ratios of $pval_{without cutting}/pval_{with cutting}$')
            ax.set_yscale('log')
            ax.set_xscale('log')
            i += 1
        pp.savefig()
        
    fig = plt.figure(figsize=(23.38, 16.54), dpi=300)
    with PdfPages('pvals_GREcut_vs_bothcut.pdf') as pp:
        i = 1
        for key in ['all','out','in','in_ex']:
            ax = fig.add_subplot(2,2,i)
            df = lists[key]
            df_filt = df.loc[df.pn <= 0.01]
            ax.scatter(df.pv_GRE_cut, df.pv_both_cut, s = 0.2,color = '#deebf7')
            ax.scatter(df_filt.pv_GRE_cut, df_filt.pv_both_cut, s = 0.2,color = '#3182bd')
            ax.set_xlabel('GRE cut')
            ax.set_ylabel('GRE & TF cut')
            ax.set_title('$pval_{without cutting}/pval_{with cutting}$')
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylim([min(lists[key].pv_both_cut), 1])
            ax.set_xlim([min(lists[key].pv_GRE_cut), 1])
            i += 1
        pp.savefig()       




#the best n% overall matches

import matplotlib_venn as venn
overlap_sets = {}
for p in [0.0001,0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
    toplists = {}
    overlap_sets[p] = {}
    for key in lists.keys():
        lists[key].sort_values(by = 'pv', inplace = True)
        toplists[key] = lists[key].iloc[0:ceil(p*len(lists[key]))][['GRE_ID', 'TF_ID', 'pv']]
        toplists[key] = list(zip(zip(toplists[key].GRE_ID, toplists[key].TF_ID), toplists[key].pv))

    pp = PdfPages('overlap_of_matches_motif_vs_position_upper_{:.2f}percent.pdf'.format(p*100))
    f = plt.figure(figsize=(23.38, 16.54), dpi=300)
    f.suptitle('Comparison for sets of TF-GRE mappings\n upper {:.2f} %'.format(p*100))
    big_ax = f.add_subplot(1,1,1, frameon = False)
    big_ax.set_ylim(0,6)
    big_ax.set_yticks(np.arange(0.5,6,1))
    big_ax.set_yticklabels(['positions', 'genes','all','out', 'in_ex', 'in'])
    big_ax.set_xlim(0,6)
    big_ax.set_xticks(np.arange(0.5,6,1))
    big_ax.set_xticklabels(['in', 'in_ex','out', 'all', 'genes', 'positions'])
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 20)
    i = 0
    for key1 in toplists.keys():
        if key1 == 'genes':
            key1_ = 'position: genes'
        elif key1 == 'positions':
            key1_ = 'position: position'
        else:
            key1_ = 'motifs: {}'.format(key1)
        for key2 in toplists.keys():
            if key2 == 'genes':
                key2_ = 'position: genes'
            elif key2 == 'positions':
                key2_ = 'position: position'
            else:
                key2_ = 'motifs: {}'.format(key2)
            i += 1
            if key1 == key2:
                continue
            else:
                ax = f.add_subplot(6,6,i)
                set1 = set([x[0] for x in toplists[key1]])
                set2 = set([x[0] for x in toplists[key2]])
                overlap = set1.intersection(set2)
                overlap_sets[p][(key1,key2)] = overlap
                v = venn.venn2([set([x[0] for x in toplists[key1]]), set([x[0] for x in toplists[key2]])], set_labels=(key1_, key2_))  
    plt.subplots_adjust(left  = 0.1, right = 0.9, bottom = 0.1, top = 0.9, wspace = 1, hspace = 0.2)
    pp.savefig()
    pp.close()




sublists_GRE = {} # a dictionary: for all GREs a list of the matched TFs
sublists_TF = {} # a dictionary: for all TFs a list of the matched GREs
TF_IDs = set(lists['all'].TF_ID)
GRE_IDs = set(lists['all'].GRE_ID)
progbar = pb.ProgressBar()
for key in progbar(lists.keys()):
    print(key)
    sublists_GRE[key] = {} 
    sublists_TF[key] = {}
    for GI in GRE_IDs:
        sublists_GRE[key][GI] = lists[key].loc[lists[key].GRE_ID == GI]
        sublists_GRE[key][GI].sort_values(by = 'pv', inplace = True)        
    for TF in TF_IDs:
        print(TF)
        sublists_TF[key][TF] = lists[key].loc[lists[key].TF_ID == TF]
        sublists_TF[key][TF].sort_values(by = 'pv', inplace = True)



#for every TF [GRE], find the n best matching GREs [TFs] for each of the mapping approaches
top_n_TFs = {}
top_n_GREs = {}
ppp = pb.ProgressBar()
for typ in ppp(sublists_GRE.keys()):
    top_n_TFs[typ] = {}
    top_n_GREs[typ] = {}
    for i in range(1,11):
        top_n_TFs[typ][i] = {}
        top_n_GREs[typ][i] = {}
        for G in GRE_IDs:
            top_n_TFs[typ][i][G] = set(sublists_GRE[typ][G].iloc[0:i].TF_ID)
        for T in TF_IDs:
            top_n_GREs[typ][i][T] = set(sublists_TF[typ][T].iloc[0:i].GRE_ID)


#for every TF [GRE], find all GREs [TFs] which are among the n best matches in any of the mappings, then assign a unique value to each of the GREs [TFs] indicating which mappings support that match
multimap_GRE = {}
multimap_TF = {}
ppp = pb.ProgressBar()
for i in ppp(range(1,11)):
    multimap_GRE[i] = {}
    multimap_TF[i] = {}
    the_TFs = {}
    the_GREs = {}
    for G in GRE_IDs:
        multimap_GRE[i][G] = {}
        the_TFs = set()
        for typ in top_n_TFs.keys():
            TFs = top_n_TFs[typ][i][G] 
            the_TFs = the_TFs.union(TFs) 
        for T in the_TFs:
            index = 0
            if T in top_n_TFs['in_ex'][i][G]:
                index += 1
            if T in top_n_TFs['in'][i][G]:
                index += 2
            if T in top_n_TFs['out'][i][G]:
                index += 4
            if T in top_n_TFs['all'][i][G]:
                index += 8
            if T in top_n_TFs['genes'][i][G]:
                index += 16
            if T in top_n_TFs['positions'][i][G]:
                index += 32
            multimap_GRE[i][G][T] = index
    for T in TF_IDs:
        multimap_TF[i][T] = {}
        the_GREs = set()
        for typ in top_n_GREs.keys():
            GREs = top_n_GREs[typ][i][T] 
            the_GREs = the_GREs.union(GREs) 
        for G in the_GREs:
            index = 0
            if G in top_n_GREs['in_ex'][i][T]:
                index += 1
            if G in top_n_GREs['in'][i][T]:
                index += 2
            if G in top_n_GREs['out'][i][T]:
                index += 4
            if G in top_n_GREs['all'][i][T]:
                index += 8
            if G in top_n_GREs['genes'][i][T]:
                index += 16
            if G in top_n_GREs['positions'][i][T]:
                index += 32
            multimap_TF[i][T][G] = index

#find high confidence matches on three levels
for Tkey in multimap_TF[1].keys():  
    for Gkey, i in multimap_TF[1][Tkey].items():
        if i == 63:
            print( '{} and {} are a perfect match based on best matches for TFs!!'.format(Tkey, Gkey))
        if i in [62,61,59,55,47,31]:
            print( '{} and {} are almost a perfect match based on best matches for TFs!!'.format(Tkey, Gkey))
        if i in [60,58,54,46,30,57,53, 45,29,51,43,27,39,23,15]:
            print( '{} and {} are a decent match based on best matches for TFs!'.format(Tkey, Gkey))
for Gkey in multimap_GRE[1].keys():  
    for Tkey, i in multimap_GRE[1][Gkey].items():
        if i == 63:
            print( '{} and {} are a perfect match based on best matches for GREs!'.format(Tkey, Gkey))
        if i in [62,61,59,55,47,31]:
            print( '{} and {} are almost a perfect match based on best matches for GREs!'.format(Tkey, Gkey))
        if i in [60,58,54,46,30,57,53, 45,29,51,43,27,39,23,15]:
            print( '{} and {} are a decent match based on best matches for GREs!'.format(Tkey, Gkey))
        



#extract high confidence matches from gene/position mapping
            
df = pd.merge(ml_g, ml_p, on=['GRE_ID', 'TF_ID'], suffixes = ['_genes', '_positions'])
alpha = 0.05
bf = alpha/len(df) 
based_on_gene_set= set([(x.GRE_ID, x.TF_ID) for _, x in df.loc[df.pv_genes <bf].iterrows()])
based_on_position_set = set([(x.GRE_ID, x.TF_ID) for _, x in df.loc[df.pv_positions <bf].iterrows()])
based_on_both_set = based_on_gene_set.intersection(based_on_position_set)  


#create tomtom output to visually evaluate motif similarity for high confidence position matches
from tomtom_for_candidates import tomtom_candidates
for ps in ['all', 'out', 'in', 'in_and_expressed']: #matches significant on gene and position level
    for p in ['pal', 'nonpal']:
        for match in based_on_both_set:
            tomtom_candidates(GRE_ID = match[0], peak_set = ps, TF_ID = match[1], TF_motif_shape = p, e_thresh = 1, dist_method = 'pearson', min_overlap = 5)



#plot pvalues of gene mapping vs position mapping
pp = PdfPages('pval_genes_vs_pval_position_position_mapping.pdf')
df = pd.merge(ml_g, ml_p, on=['GRE_ID', 'TF_ID'], suffixes = ['_genes', '_positions'])
alpha = 0.05
bf = alpha/len(df) 
fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
plt.plot(df.pv_genes, df.pv_positions, 'k.',markersize = 1, label = "all pvalue pairs")
plt.xlabel('pvalue based on genes')
plt.ylabel('pvalue based on position')
plt.vlines(alpha, ymin = min(df.pv_positions), ymax = 1.5, colors = "r", linestyles = 'dashed', linewidth = 1, label = '$p = \\alpha = %s$'%alpha)
plt.vlines(bf, ymin = min(df.pv_positions), ymax = 1.5, colors = "r", linestyles = 'dotted', linewidth = 1, label = '$p = \\alpha_{BF}$ = %3.2e'%bf)
plt.hlines([alpha,bf], xmin = min(df.pv_genes), xmax = 1.5, colors = "r", linestyles = ['dashed', 'dotted'], linewidth = [1,1])
plt.yscale('log')
plt.xscale('log')
plt.xlim([min(df.pv_genes)/10, 5])
plt.ylim([min(df.pv_positions)/10, 5])
n_alpha_genes = len(df.loc[df.pv_genes < alpha])
n_alpha_positions = len(df.loc[df.pv_positions < alpha])
n_alpha_both = len(df.loc[df.pv_genes < alpha].loc[df.pv_positions < alpha])
n_bf_genes = len(df.loc[df.pv_genes < bf])
n_bf_positions = len(df.loc[df.pv_positions < bf])
n_bf_both = len(df.loc[df.pv_genes < bf].loc[df.pv_positions < bf])
plt.text(x = 0.1, y = 0.5, s = ' {} matches below gene threshold, \n {} matches below position threshold, \n {} matches below both thresholds'.format(n_bf_genes, n_bf_positions, n_bf_both), transform = plt.gca().transAxes)
plt.title('pvalues for matches based on gene mapping and position mapping')
plt.legend(loc = 3)
pp.savefig()
pp.close()
plt.close(fig)

#the sets of significant matches:
based_on_gene_set= set([(x.GRE_ID, x.TF_ID) for _, x in df.loc[df.pv_genes <bf].iterrows()])
based_on_position_set = set([(x.GRE_ID, x.TF_ID) for _, x in df.loc[df.pv_positions <bf].iterrows()])
based_on_both_set = based_on_gene_set.intersection(based_on_position_set)

cw = csv.writer(open('matches_based_on_both_mappings.csv','w'))
for match in based_on_both_set:
    cw.writerow(match)
cw.close()
###---------------------------------------------------------------------------------------------###            

##create a mapping of TF to genes via TF -- genomic position -- gene promoter region
#
#
##dictionary of 3996 genes with values being the region - 70 + 150 bp from TSS (promoter region)
#gene_dic = pickle.load(open('Dictionaries/gene_dic.p', 'br'))
##dictionary of 156 TFs with values being a list of positions (+- 12bp around the binding position)
#chip_dic = pickle.load(open('Dictionaries/chip_dic.p', 'br'))
#
#chip_gene_promoters = {}
#pbar = pb.ProgressBar()
#for TF, sites in pbar(chip_dic.items()):
#    chip_gene_promoters[TF] = []
#    for site in sites:  
#        central = mean(site)
#        for gene, promoter in gene_dic.items():
#            if (promoter[0] <= central <= promoter[1]):
#                chip_gene_promoters[TF].append(gene)              #consider a TF-gene match if the central binding position of TF lies inside the promoter region of gene
#                
#
#pickle.dump(chip_gene_promoters, open('Dictionaries/chip_gene_promoters.p', 'bw'))
#                
#chip_gene_promoters = pickle.load(open('Dictionaries/chip_gene_promoters.p', 'br'))
#

###---------------------------------------------------------------------------------------------###
          

#create a mapping of TF to genes via TF -- genomic position -- gene region


##extract dictionary of 3996 genes with values being genomic region from csv file
#os.chdir('/Users/rkoch/Documents/Data_and_Scripts')
#with open('./Data/gene_data.csv', 'r') as gd:
#    gene_data = gd.readlines()
#gene_coordinates_dic = {}
#for row in gene_data[1:]:
#    row = row.split(',')
#    gene_coordinates_dic[row[0]] = [int(x) for x in row[1:3]]
#pickle.dump(gene_coordinates_dic, open('Dictionaries/gene_coordinates_dic.p', 'bw'))
#gene_coordinates_dic = pickle.load(open('Dictionaries/gene_coordinates_dic.p', 'br'))
##dictionary of 156 TFs with values being a list of positions (+- 12bp around the binding position)
#chip_dic = pickle.load(open('Dictionaries/chip_dic.p', 'br'))
#
#chip_genes = {}
#pbar = pb.ProgressBar()
#for TF, sites in pbar(chip_dic.items()):
#    chip_genes[TF] = []
#    for site in sites:  
#        central = mean(site)
#        for gene, pos in gene_coordinates_dic.items():
#            if (pos[0] <= central <= pos[1]):
#                chip_genes[TF].append(gene)              #consider a TF-gene match if the central binding position of TF lies inside the genomic region of gene
#                
#
#pickle.dump(chip_genes, open('Dictionaries/chip_genes.p', 'bw'))
#                
#chip_genes = pickle.load(open('Dictionaries/chip_genes.p', 'br'))


###---------------------------------------------------------------------------------------------###
#
##compare the sources for TF-gene mapping:
#from matplotlib import pyplot as plt
#import matplotlib_venn as venn
#pcd = set(pub_chip_dic.keys())
#cd = set(chip_dict.keys())
#cg = set(chip_gene_promoters.keys())
#
#v = venn.venn3_unweighted([pcd, cd, cg], set_labels=('pub_chip_dic', "chip_dict", 'chip_gene_promoters'))
#difference_TF_sets = pcd.difference(cd)
#
#
#for key in cd:
#    targets_cd = set(chip_dict[key])
#    size_cd = len(targets_cd)
#    targets_pcd = set(pub_chip_dic[key])
#    size_pcd = len(targets_pcd)
#    targets_cg = set(chip_gene_promoters[key])
#    size_cg = len(targets_cg)
#    plt.figure()
#    venn.venn2([targets_pcd, targets_cg], set_labels =('pub_chip_dic', 'chip_gene_promoters'))
#    plt.title('targets for TF %s'%key)
#    plt.figure()
#    plt.bar(['pub_chip_dic', 'chip_gene_promoters', 'chip_dict'], [size_pcd, size_cg, size_cd])
#
#
#total_targets_pcd = sum([len(x) for x in pub_chip_dic.values()])
#total_targets_cd = sum([len(x) for x in chip_dict.values()])
#total_targets_cg = sum([len(x) for x in chip_gene_promoters.values()])

###---------------------------------------------------------------------------------------------### 

#compare logos for TFs based on different peak sets
import subprocess
os.chdir('/Users/rkoch/Documents/Data_and_Scripts/MEME_output')
no_of_genes = {}
types = ['in', 'out', 'in_and_expressed', 'all']
TFs = [x for x in os.listdir('in') if x.startswith('Rv')]
shapes = ['pal','nonpal']
TFs = [x.split("_")[0] for x in TFs]
emptyfiles = {'in':0, 'out':0, 'in_and_expressed':0, 'all':0}
for TF in TFs:
    pp = PdfPages('compare/%s.pdf'%TF)
    i = 0
    fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
    no_of_genes[TF] = {}
    for typ in types:
        no_of_genes[TF][typ] = {}
        for shape in shapes:
            no_of_genes[TF][typ][shape] = 0
            memefile = ('%s/%s_%s/meme.txt'%(typ, TF, shape))
        #extract the number of genes used for creation of motif
            with open(memefile) as f:
                lines = f.readlines()
                if len(lines) > 0:
                    line1 = [x for x in lines if x.startswith('data:')] #line containing number of genes
                    no = line1[0].split(" ")[-1].split('\n')[0]
                    no_of_genes[TF][typ][shape] = int(no)
                else:
                    emptyfiles[typ] += 1
                outdir = ('temp_out')
                command = ['/usr/local/bin/meme2images',
                           memefile,
                           outdir]
                a = subprocess.run(command, check = True)
                print(a.returncode)
                for j in [1,2]:
                    command2 = ['/usr/local/Cellar/imagemagick/7.0.7-28/bin/convert',               #WHY would subprocess return 1 but it works in bash?
                                './temp_out/logo%s.eps'%j, 
                                './temp_out/logo%s.png'%j]
                    b = subprocess.run(command2, check = True)
                    print(b.returncode)
                    i+=1
                    ax = fig.add_subplot(4,4,i)
                    ax.imshow("./temp_out/logo%s.png"%j)
                    ax.text()
                    ax.title('%s%s%s%s'%(TF,typ,shape, j))
    pp.savefig()
###---------------------------------------------------------------------------------------------###            

#check the pvalues for recall of genes via FIMO
os.chdir('/Users/rkoch/Documents/Data_and_Scripts')
no_exp_binding_TFs_dic = pickle.load(open('./TFs_without_exp_binding.p', 'rb'))
AUPR_dic = pickle.load(open('./AUPR_dic_FIMO_vs_exp.p', 'rb'))
pval_rs_dic = pickle.load(open('./pvals_ranksum_FIMO_vs_exp.p', 'rb'))
pval_AUPR_dic = pickle.load(open('./pvals_AUPR_FIMO_vs_exp.p', 'rb'))

with PdfPages('AUPR_and_ranksum_pvalues.pdf') as pdfx:
    for sett in AUPR_dic.keys():
        plt.figure(figsize=(11.69, 8.27), dpi=100)
        ax = plt.subplot(121)
        x = plt.hist(pval_AUPR_dic[sett].values(), bins = np.arange(0,1.025,0.025))
        plt.title('Distribution of AUPR pvalues for set %s'%sett)
        plt.ylabel('frequency')
        plt.xlabel('p value')
        plt.plot([0.05,0.05], [0,max(x[0])], '--r')
        ax = plt.subplot(122)
        x = plt.hist(pval_rs_dic[sett].values(), bins = np.arange(0,1.025,0.025))
        plt.title('Distribution of ranksum pvalues for set %s'%sett)
        plt.xlabel('p value')
        plt.plot([0.05,0.05], [0,max(x[0])], '--r')
        pdfx.savefig()

###---------------------------------------------------------------------------------------------###            
#Repositories: Corem-Gene-assignments etc

#dictionary of 3996 genes with values being the region - 70 + 150 bp from TSS (promoter region)
gene_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/gene_dic.p', 'br'))


#dictionary of 156 TFs with values being a list of positions (+- 12bp around the binding position)    
chip_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_dic.p', 'br'))


#dictionary of 143 TFs with values being a list of genes bound by them
chip_dict = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_dict.p', 'br'))
#missing TFs: Rv3334, Rv1963c, Rv0821c, Rv1152, Rv2703, Rv2160c, Rv0020c, Rv1626, Rv2258c, Rv1151c, Rv0377, Rv1909c, Rv3060c


#dictionary of 156 TFs with values being a list of genes (from gene_dic * chip_dic)
chip_gene_promoters = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_gene_promoters.p', 'br'))
#alternatively (whole +-12bp window must be in the promoter region) -> a little less links
#chip_gene = pickle.load(open('Dictionaries/chip_gene.p', 'br'))


#dictionary of 209 TFs with entries being dictionaries of all genes as keys and foldchange expression change as values
tf_oe_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/tf_oe_dic.p', 'br'))


#dictionary of 156 TFs with entries being lists of genes bound by them. Probably from public resources?
pub_chip_dic = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/pub_chip_dic.p', 'br'))


###---------------------------------------------------------------------------------------------###            
#build sets of TFs etc

all_genes_set = set([x for x,z in gene_dic.items() if z])

TFs_with_peaks_set = set([x for x,z in chip_dic.items() if z])

TFs_with_gene_hits_pcd_set = set([x for x,z in pub_chip_dic.items() if z])

TFs_with_gene_hits_cd_set = set([x for x,z in chip_dict.items() if z])

TFs_overexpressed_set = set([x for x,z in tf_oe_dic.items() if z])




###---------------------------------------------------------------------------------------------###            

#Robin creates his own repositories from scratch from the ChipSeqData and other sources

#base repository: publication ChipSeq data Minch et al
big_df = pd.read_csv(open('/Users/rkoch/Desktop/literature/minch/ncomms6829-s2.csv', 'r'))
all_TFs = set(big_df.Regulator)

#a dictionary containing for every TF the peak centers where it bound in the ChipSeq data
TF_peaks = {}
for TF in all_TFs:
    peaks = [int(x) for x in big_df.loc[big_df.Regulator == TF]['Peak Center (Ccenter)']]
    TF_peaks[TF] = list(peaks)
pickle.dump(TF_peaks, open('TF_peaks.p', 'wb'))

#extract dictionary of 3996 genes with values being genomic region from csv file
os.chdir('/Users/rkoch/Documents/Data_and_Scripts')
with open('./Data/gene_data.csv', 'r') as gd:
    gene_data = gd.readlines()
gene_coordinates_dic = {}
for row in gene_data[1:]:
    row = row.split(',')
    gene_coordinates_dic[row[0]] = [int(x) for x in row[1:3]]
pickle.dump(gene_coordinates_dic, open('Dictionaries/gene_coordinates_dic.p', 'bw'))


#a dictionary of TFs with values being the genes where TF bound the gene sequence
TF_genebodies_dic = {}
TF_genes_dic = {}
probar = pb.ProgressBar()
for TF in probar(all_TFs):
    TF_genes_dic[TF] = []
    for peak in TF_peaks[TF]:
        for gene, positions in gene_coordinates_dic.items():
            if positions[0] < peak < positions[1]:
                TF_genebodies_dic[TF].append(gene) 


#I would not trust this dic, since who knows whether position 1 of a gene is also their TSS?                
#gene_promoters_dic = {}
#for gene, positions in gene_coordinates_dic.items():
#    gene_promoters_dic[gene] = [positions[0]-70, positions[0]+150]



#a dictionary of TFs with values being the genes where TF bound the promoter
pbb = pb.ProgressBar()
TF_promoter_bound_genes = {}
for TF, peaks in pbb(TF_peaks.items()):
    TF_promoter_bound_genes[TF] = []
    for peak in peaks:
        for gene, boundaries in gene_dic.items():
            if boundaries[0] < peak < boundaries[1]:
                TF_promoter_bound_genes[TF].append(gene)
pickle.dump(TF_promoter_bound_genes, open('TF_promoter_bound_genes.p', 'wb'))

#create a dict of all genes diff expressed upon OE of TF  
TF_gene_diffex = pd.read_csv(open('/Users/rkoch/Documents/Data_and_Scripts/Data/TF_gene_diffex.csv', 'r'))      #derived from ncomms6829-s4 Minch et al   
TF_diffex_dic = {}
TF_diffex_dic = {k:[] for k in set(TF_gene_diffex.Regulator)}
for _, row in TF_gene_diffex.iterrows():
    if row.diffex in ['IND', 'REP']:
        TF_diffex_dic[row.Regulator].append(row.gene)
pickle.dump(TF_diffex_dic, open('TF_diffex_dic.p', 'wb'))
        
        
        
TF_diffex_bound = {k:[] for k in TF_promoter_bound_genes.keys()}
for TF, genes in TF_diffex_dic.items():
    for gene in genes:
        if gene in TF_promoter_bound_genes[TF]:
            TF_diffex_bound[TF].append(gene)
pickle.dump(TF_diffex_bound, open('TF_diffex_bound.p', 'wb'))


    

#create corem-condition dict
    
corem_cond_df = pd.read_csv(open("/Users/rkoch/Documents/Data_and_Scripts/Data/corem_cond.csv","r"))

all_corems = set(corem_cond_df.COREM)

corem_cond_dic = {}
for corem in all_corems:
    cdf = corem_cond_df.loc[corem_cond_df.COREM == corem]
    corem_cond_dic[corem] = list(cdf["EGRIN2.block"])

pickle.dump(corem_cond_dic, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/corem_cond_dic.p', 'wb'))
    