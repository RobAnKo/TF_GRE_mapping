#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 14:46:37 2018

@author: rkoch

precision-recall for FIMO vs experimental TF binding results

"""
import pickle
import os
import subprocess
import myfuncs as mf
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import progressbar as pb
from matplotlib.backends.backend_pdf import PdfPages
nan = float('NaN')
    
def run_fimo(TF_name, peakset = 'in', palindromic = False, target = 'promoter', thresh = 0.05, output = 'genes'):
    """
    a wrapper to run the FIMO tool and collect the data, either genes or positions
    
    
    """   
    if not TF_name.startswith('Rv'):
        print("Error: TF_name must be a valid ID, eg 'Rv0022c'")
        return
    if peakset not in ['in', 'out', 'in_and_expressed', 'all']:
        print("Error: peakset has to be one of 'in', 'out', 'in_and_expressed' and 'all'")
        return
    if type(palindromic) != bool:
        print('Error: palindromic has to be boolean')
        return
    if target not in ['promoters', 'corem_promoters', 'genes']:
        print("Error: target has to be either 'promoters' or 'corem_promoters' or 'genes'")
        return
    if not mf.is_number(thresh):
        print("Error: threshold has to be a number between 0 and 1")
        return
    if not output in ['genes', 'positions']:
        print("Error: output has to be either 'genes' or 'positions'")
        return
    else:
        thresh = float(thresh)
        
    
    if target == 'promoters':
        bgfile = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters_bgmodel.txt'
        database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters.fna'
    elif target == 'corem_promoters':
        bgfile = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/corem_gene_promoters_bgmodel.txt'
        database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/corem_gene_promoters.fna'
    else:
        bgfile = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/genes_bgmodel.txt'
        database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/genes.fna'

    if palindromic:
        p = '_pal'
    else:
        p = '_nonpal'

    file = '/Users/rkoch/Documents/Data_and_Scripts/MEME_output/%s/%s%s/meme.txt'%(peakset, TF_name, p)  

    command = ['/Users/rkoch/Programs/meme/bin/fimo', 
               '--bgfile', bgfile,
               '--oc', 'temp_out', 
               '--norc',
               '--qv-thresh', 
               '--motif', str(1),
               '--thresh', str(thresh),
               '--max-stored-scores', str(1e5),#reduce hits to 100000 to speed up analysis
               '--verbosity', str(1),
               '--parse-genomic-coord',
               file, 
               database]

    subprocess.run(command)
    with open('temp_out/FIMO.txt', 'r') as f:
        result = f.readlines()
    result.pop(0)
    out = [nan]*len(result)
    i = 0
    
    if output == 'genes':
        for row in result:
            row = row.split('\t')
            name = row[2]
            out[i] = name
            i+=1
        out = mf.uniquify(out)
        print('%s genes were found'%len(out))
    else:
        for row in result:
            row = row.split('\t')
            pos = [int(x) for x in row[3:5]]
            out[i] = pos
            i+=1
        print('%s positions were found'%len(out))         
    return out       



def overlap_positions(boundaries1, list_boundaries2, t):
    """
    This function takes as input two positional intervals as pairs of boundary positions 
    and decides whether the intervals overlap in any way
    
    """
    b1low, b1up = min(boundaries1), max(boundaries1)
    for boundaries2 in list_boundaries2:
        b2low, b2up = min(boundaries2), max(boundaries2)
        if not (b2low > b1up or b1low > b2up):                        #cases of no overlap
            return True, list_boundaries2.index(boundaries2)
    return False, None

def overlap_genes(gene1, list_genes2):
    if list_genes2.count(gene1):
        return True, list_genes2.index(gene1)
    return False, None


def polygon_area(x,y):
    correction = x[-1] * y[0] - y[-1]* x[0]
    main_area = np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:])
    return 0.5*np.abs(main_area + correction)



def calculate_AUPR(predictions, data, AUPR_only = False, typ = 'genes'):
    """
    This function calculates the AUPR as well as optionally precision and recall values for a sorted 
    set of predicted items compared to a validated set of items.
    
    input:
        - predictions: a list of length n of predicted items sorted from high to low confidence predicitions
        - data: a list of validated items
        - AUPR_only: boolean. if True, only one float value for AUPR is returned, else precision and recall are returned as lists of length n 
    output:
        a float for AUPR, optionally precision and recall values along the predictions
    """
    curve = [nan]*len(predictions)
    N = 0 #number of predictions
    P = len(data)
    TP = 0
    #progner = pb.ProgressBar()
    for item in predictions:#progner(predictions):
        N += 1
        if typ == 'genes':
            ol = overlap_genes(item, data)
        if typ == 'positions':
            ol = overlap_positions(item, data, t = 5)
        if ol[0]:
            TP += 1
            data.pop(ol[1])
            prec = TP/N
            reca = TP/P
            curve[N-1] = [prec,reca]
    curve = [x for x in curve if x is not nan]
    if curve:
        precision = [x[0] for x in curve]
        recall = [x[1] for x in curve]
        precision.extend([0,0,precision[0]])
        recall.extend([recall[-1], 0, 0])
        area = polygon_area(precision, recall)
    else:
        area = None
        precision = [None]
        recall = [None]
    if AUPR_only:
        return area
    else:
        return area, precision, recall


os.chdir('/Users/rkoch/Documents/prec_rec_output/')
exp_peaks = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_dic.p', 'rb'))

##import positions of all genes as well as all genes bound by each TF to produce a list of all positions of genes bound by TF
#with open('/Users/rkoch/Documents/Data_and_Scripts/Data/gene_data.csv', 'r') as f:
#    lines = f.readlines() 
#gene_positions = {}                                                          #gene: positions
#for line in lines[1:]:
#    line = line.split(',')
#    gene = line[0]
#    start = int(line[1])
#    end = int(line[2])
#    gene_positions[gene] = [start, end]                                    
exp_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/chip_dict.p', 'rb'))           #TF: genes bound
#exp_genes_positions = {}                                                    #TF: gene positions
#for key in exp_genes.keys():
#    exp_genes_positions[key] = []
#    for gene in exp_genes[key]:
#        exp_genes_positions[key].append(gene_positions[gene])

for typ in ['genes', 'positions']:
#level1: type of experimental data: genes vs positions
    AUPR_dic = {}
    no_exp_binding_TFs_dic = {}
    #pval_rs_dic = {}
    pval_AUPR_dic = {}
    for peak_set in ['all','out','in','in_and_expressed']:    
    #level2: set of peaks used to create motifs  
        no_exp_binding_TFs_dic[peak_set] = {}
        AUPR_dic[peak_set] = {}
        pval_AUPR_dic[peak_set] = {}
        for p in (True, False):
        #level3: palindromic vs non-palindromic motifs
            if p:
                pbit = '_pal'
            else:
                pbit = '_nonpal'
            ##produce precision recall curves for overlap of genes experimentally shown to be bound by TF vs genes containing FIMO hits with the derived TF motif
            TFs = list(exp_genes.keys())
            AUPRs = {}
            #pvals_rs = {}
            pvals_AUPR = {}
            no_exp_binding_TFs = []
            with PdfPages('PRCs_%s_%s%s.pdf'%(typ, peak_set, pbit)) as pdf:  
                print('working on set %s%s using %s'%(peak_set, pbit, typ))
                progger = pb.ProgressBar()
                for TF in progger(TFs):
                #level4: single TFs
                    print(TF)
                    if typ == 'genes':
                        data = exp_genes[TF]
                    elif typ == 'positions':
                        data = exp_peaks[TF]
                    if not data:
                        print("no PRC could be generated for TF %s because no experimental validated binding was found"%TF)
                        no_exp_binding_TFs.append(TF)
                        AUPRs.append(None)
                        continue
                               
                    #produce FIMO hits
                    #if peak_set in ['in', 'in_and_expressed']:
                    predictions = run_fimo(TF, peakset = peak_set, palindromic = p, target = 'promoters', thresh = 1, output = typ)
                    #elif peak_set in ['all', 'out']:
                    #    predictions = run_fimo(TF, peakset = peak_set, palindromic = p, target = 'genes', thresh = 1, output = typ)
                    
                        
                    AUPR, precision, recall = calculate_AUPR(predictions, data[:], typ = typ) #copy list to avoid mutation of original list since it is needed later
                    AUPRs[TF] = AUPR
                    #plot
                    plt.figure(figsize=(11.69, 8.27), dpi=100)
                    if AUPR == None:
                        plt.text(x = 0.05, y = 0.5, s='no PRC could be generated for TF %s \n because no experimental validated binding was found'%TF)
                        continue
                    else:
                        plt.xlim((0,1))
                        plt.ylim((0,1))
                        plt.xlabel('recall')
                        plt.ylabel('precision')
                        plt.title('precision-recall curve for recovering TF-bound {} using FIMO with {}{} motifs \nfor TF {} AUPR: {:.3f}'.format(typ, peak_set, pbit, TF, AUPR))
                        ratio = len(data)/len(predictions)
                        plt.plot(recall, precision, '-k', label = "PRC")
                        plt.plot([0,1],[ratio,ratio],'--r', label = "random predictor")
                        plt.legend()
                        pdf.savefig()
                        
                    
                    #found_hits = [x for x in exp_hits if x in  predictions]
                    #calculate rank sum of predictions
                    #ranksum = sum([predictions.index(x) for x in found_hits])
                    #calculate a null distribution of ranksums to get pvalue
                    copy = predictions[:]
                    #ranksum_distribution = [nan]*10000
                    AUPR_distribution = [nan]*10000
                    print('creating null distributions for AUPRs and rank sums')
                    prog = pb.ProgressBar()
                    for i in prog(range(10000)):
                        #print(i)
                        np.random.shuffle(copy)
                        #ranksum_distribution[i] = sum([copy.index(x) for x in found_hits])
                        AUPR_distribution[i] = calculate_AUPR(copy, data[:], AUPR_only = True, typ = typ)
                    #pval_rs = len([x for x in ranksum_distribution if x < ranksum])/10000
                    pval_AUPR = len([x for x in AUPR_distribution if x > AUPR])/10000
                    #pvals_rs[TF] = pval_rs
                    pvals_AUPR[TF] = pval_AUPR
                    
              
                    
            no_exp_binding_TFs_dic[peak_set][pbit] = no_exp_binding_TFs
            AUPR_dic[peak_set][pbit] = AUPRs    
            #pval_rs_dic[(peak_set, pbit)] = pvals_rs
            pval_AUPR_dic[peak_set][pbit] = pvals_AUPR
    pickle.dump(no_exp_binding_TFs_dic, open('./TFs_without_exp_binding_%s.p'%typ, 'wb'))
    pickle.dump(AUPR_dic, open('./AUPR_dic_FIMO_vs_exp_%s.p'%typ, 'wb'))
        #pickle.dump(pval_rs_dic, open('./pvals_ranksum_FIMO_vs_exp.p', 'wb'))
    pickle.dump(pval_AUPR_dic, open('./pvals_AUPR_FIMO_vs_exp_%s.p'%typ, 'wb'))
