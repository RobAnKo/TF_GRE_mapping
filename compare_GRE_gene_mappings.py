#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 14:00:09 2018

@author: rkoch
"""

import pandas as pd
import pickle 
import os
import matplotlib_venn as venn
import subprocess

####----------------------------0. find a list of all GREs--------------------------------------####

match_list = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/combined_match_list_all.p', 'rb'))

unique_GREs = set([x[0] for x in match_list])
####-----------I. create GRE-gene connection via GRE-corem--corem-gene datasets-----------------####
GRE_corems = pd.read_csv('/Users/rkoch/Documents/Data_and_Scripts/Data/gre_corem.csv')
GRE_corems.gre = GRE_corems.gre -1

corems_genes = pickle.load(open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/corems_genes.p', 'rb'))



unique_corems = set(GRE_corems.corem)

GRE_cor_dic = {}
for GRE in unique_GREs:
    GRE_cor_dic[GRE] = GRE_corems.loc[GRE_corems.gre == GRE].corem.values
    
GRE_cor_genes = {}

for GRE in unique_GREs:
    GRE_cor_genes[GRE] = []
    corems = GRE_cor_dic[GRE]
    for corem in corems:
        GRE_cor_genes[GRE].extend(corems_genes[corem])
        
for GRE, genes in GRE_cor_genes.items():
    GRE_cor_genes[GRE] = set(genes)
  
pickle.dump(GRE_cor_genes, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/GRE_cor_genes.p', 'wb'))
 



####-------------------II. load GRE-gene connections from files--------------------------------####
    
dir_ = '/Users/rkoch/Documents/Data_and_Scripts/Data/gene_all_gres/'
GRE_genes = {}
for GRE in unique_GREs:
    GRE_genes[GRE] = set()
    ID = GRE+1
    file = dir_ + 'GRE{}.csv'.format(ID)
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = [x.split('"') for x in f.readlines()]
            genes = [x[1] for x in lines[1:] if len(x)>1]
            GRE_genes[GRE] = set(genes)
pickle.dump(GRE_genes, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/GRE_genes.p','wb'))




####----------------III. run FIMO for every GRE against promoters and find hits------------------####

GRE_FIMO_genes = {}
bg_file = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters_bgmodel.txt'
thresh = 1
database = '/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/gene_promoters.fna'
snap = pb.ProgressBar()
for GRE in snap(unique_GREs):
    GRE_FIMO_genes[GRE] = set()
    GRE_file = '/Users/rkoch/Documents/Data_and_Scripts/motif_clusters_24/{:04.0f}_memeOut.txt'.format(GRE)   
    no_hits = max((len(GRE_genes[GRE]), len(GRE_cor_genes[GRE])))
    GRE_command = ['/Users/rkoch/Programs/meme/bin/fimo', 
                   '--bgfile', bg_file , 
                   '--oc', 'temp_out', 
                   '--parse-genomic-coord',
                   '--qv-thresh',
                   '--thresh', '{:.10f}'.format(thresh),
                   '--max-stored-scores', str(no_hits),
                   '--verbosity', str(1),
                   GRE_file, 
                   database]
    sp = subprocess.run(GRE_command)
    
    with open('temp_out/fimo.txt') as GRE_out:
        GRE_hits = GRE_out.readlines()[1:]
        if GRE_hits:
            GRE_hits = [x.split('\t') for x in GRE_hits]
            GRE_hits = set([x[2] for x in GRE_hits])
            GRE_FIMO_genes[GRE] = GRE_hits


pickle.dump(GRE_FIMO_genes, open('/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/GRE_FIMO_genes.p','wb'))




####------------------------IV. comparison of the three mappings--------------------------------####



for GRE in unique_GREs:
    f = plt.figure()
    plt.title('{}'.format(GRE))
    venn.venn3([GRE_cor_genes[GRE], GRE_genes[GRE], GRE_FIMO_genes[GRE]], set_labels = ['GRE_cor_genes','GRE_genes', "GRE_FIMO_genes"])
    #venn.venn2([GRE_cor_genes[GRE], GRE_genes[GRE]], set_labels = ['GRE_cor_genes','GRE_genes'])