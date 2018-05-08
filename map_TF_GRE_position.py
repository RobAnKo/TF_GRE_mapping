#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 12:01:43 2018

@author: rkoch


This script generates match lists for all four possible resolutions of positional mapping:
    
    -genes
    
    -hits
    
    -nucleotides
    
    -peak centers

and saves them in the directory specified below as <out_dir>.

Prerequisites are:

    -availibility of the FIMO binary (<fimo_path>)
    
    -the set of GRE motifs stored in one folder (<GRE_directory>)
    
    - a dictionary associating TFs with the genes they bind in the promoter region
    
    - a dictionary associating TFs with the positions where they bind in TFOE experiments


"""
import subprocess
from Bio import SeqIO
import scipy.stats as ss
import pickle
import progressbar as pb
import numpy as np



genome_directory = "/Users/rkoch/Documents/Data_and_Scripts/MTB_genome_data_12/"
GRE_directory = '/Users/rkoch/Documents/Data_and_Scripts/motif_clusters_24/'
fimo_path = '/Users/rkoch/Programs/meme/bin/fimo'
data_directory = '/Users/rkoch/Documents/Data_and_Scripts/Dictionaries/'
out_dir = '/Users/rkoch/Documents/Data_and_Scripts/out/match_lists/'

TF_pos = pickle.load(open('{}chip_dic.p'.format(data_directory), 'br'))                       #repository for positions bound by TFs  
TF_promoter = pickle.load(open('{}chip_gene_promoters.p'.format(data_directory), 'br')) #repository for genes bound by TFs in their promoter region
TF_IDs = list(TF_promoter.keys())
GRE_IDs = pickle.load(open('{}GRE_IDs.p'.format(data_directory), 'rb'))
nan =np.float('NaN')    



def overlap_positions(GRE_hits, TF_hits):
    out = 0
    for Ghit in GRE_hits:
        for Thit in TF_hits:
            if Ghit[1] < Thit[0]:
                break
            if (Ghit[0] in range(Thit[0], Thit[1])) or (Ghit[0] in range(Thit[0], Thit[1])):
                out += 1
    return out
    
    
    
          
###---------------------------------main-function------------------------------------------------###     
        
def position_map_TF_GRE(GRE_ID,
                        TF_ID, 
                        thresh = 0.05, 
                        target = 'genome', 
                        level = 'genes', 
                        position_res = 'nucleotides'):
    
    """
    This function takes as input a GRE ID and a TF ID and optional parameters of analysis.
    It uses the FIMO algorithm to find significant matches of GRE motifs to  genome or promoter regions 
    and yields a significance score indicating the similarity of GRE matches and TFBS on the level of either bound genes, 
    overlapping matches or nucleotides bound by both TF and GRE 
    
    input:
        
        GRE_ID - an integer ID corresponding to the IDs given in the MTB EGRIN2 model
        
        TF_ID - a string corresponding to the TF ID
        
        thresh - a pvalue threshold for FIMO matches GRE vs background 
        
        target - either 'genome' or 'promoters' or 'corem_promoters'; sequences against which the GRE FIMO run 
        matches the motifs
              
        level - a string ['genes' or 'positions'] indicating whether matches are compared on gene or position level
        
        position_res - one of 'hits', 'nucleotides', 'peak_center' for overlap analysis on positional level or None for overlap analysis on gene level
       
    output: for the pair GRE_ID - TF_ID a list containing [GRE_ID, TF_ID, number of hits GRE, number of hits TF, number of overlapping hits, hypergeom pvalue of overlap]
        
    
    """
    
    ##check input
    if not TF_ID.startswith('Rv'):
        print("Error: TF_name must be a valid ID, eg 'Rv0022c'")
        return
    if target not in ['genome', 'promoters', 'corem_promoters']:
        print("Error: target has to be either 'genome', 'promoters' or 'corem_promoters'")
        return  
    if level not in ['genes','positions']:
        print("Error: level has to be either 'genes' or 'positions'")
        return
    if level == 'genes' and target == 'genome':
        print("Error: for overlap in bound genes, only 'promoters' and 'corem_promoters' are valid targets")
        return
    if position_res not in ['hits', 'nucleotides', 'peak_center']:
        print("Error: position_res needs to be either 'hits', 'nucleotides' or 'peak_center'" )
        return
        
    
    #define paths based on input
    if target == 'genome':
        bgfile = '{}genome_bgmodel.txt'.format(genome_directory)    #bgfiles created with fasta-get-markov from MEME suite
        database = '{}genome.fna'.format(genome_directory)
        with open(database) as f:
            n = len(SeqIO.read(f, 'fasta'))
    elif target == 'promoters':
        bgfile = '{}gene_promoters_bgmodel.txt'.format(genome_directory)
        database = '{}gene_promoters.fna'.format(genome_directory)
        with open(database) as f:
            rec = SeqIO.parse(f, 'fasta')
            rec = list(rec) 
            N = len(rec)                    #number of genes
            n = sum([len(x) for x in rec])  #number of nucleotides
    else:
        bgfile = '{}corem_gene_promoters_bgmodel.txt'.format(genome_directory)
        database = '{}corem_gene_promoters.fna'.format(genome_directory)
        with open(database) as f:
            rec = SeqIO.parse(f, 'fasta')
            rec = list(rec)
            N = len(rec)                    #number of genes
            n = sum([len(x) for x in rec])  #number of nucleotides
            
    GRE_file = '{}{:04d}_memeOut.txt'.format(GRE_directory, GRE_ID+1)
    
    

    
    
    ##FIMO run for the GRE
    GRE_command = [fimo_path, 
                   '--bgfile', bgfile, 
                   '--oc', 'temp_out', 
                   '--parse-genomic-coord',
                   '--thresh', '{:.10f}'.format(thresh),
                   '--max-stored-scores', str(1e5),
                   '--verbosity', str(1),
                   '--oc', '{}fimo_out/'.format(out_dir),
                   GRE_file, 
                   database]
    sp = subprocess.run(GRE_command)
    if sp.returncode != 0:
        print("FIMO run failed for GRE {} vs {}".format(GRE_ID, target))
        return
    
    
    ##GRE hits
    with open('{}fimo_out/fimo.txt'.format(out_dir)) as GRE_out:
        GRE_hits = GRE_out.readlines()[1:]
        GRE_hits = [x.split('\t') for x in GRE_hits]
        if GRE_hits != []:
            l = int(GRE_hits[0][4]) - int(GRE_hits[0][3]) #length of the GRE
            if level == 'genes':
                GRE_hits = set([x[2] for x in GRE_hits])
            elif level == 'positions':
                GRE_hits = [(int(x[3]),int(x[4])) for x in GRE_hits] 
                GRE_hits = set(GRE_hits)
                GRE_hits = sorted(GRE_hits, key = lambda x: x[0])
                if position_res == 'nucleotides':
                    GRE_hits = [range(x[0], x[1]+1) for x in GRE_hits]
                    GRE_hits = set([x for j in GRE_hits for x in j])
                elif position_res == 'peak_center': #only 'core' part of GRE to match against TFpeak centers instead of +-12 windows
                    GRE_hits = [range(x[0]+3, x[1]-2) for x in GRE_hits]
                    GRE_hits = set([x for j in GRE_hits for x in j])
            k_GRE = len(GRE_hits)
        else:
            k_GRE = 0
            l = 0
            if level == 'genes':
                GRE_hits = set()
            else:
                GRE_hits = []              
                
                
        
       
    ##TF hits
    if level == 'genes':
        TF_hits = set(TF_promoter[TF_ID])
        k_TF = len(TF_hits)
    else:
        TF_hits = sorted(TF_pos[TF_ID], key = lambda x: x[0])
        if position_res == 'nucleotides':
            TF_hits = [range(x[0], x[1]+1) for x in TF_hits]
            TF_hits = set([x for j in TF_hits for x in j])
        elif position_res == 'peak_center':
            TF_hits = set([int(np.mean(x)) for x in TF_hits])
        k_TF = len(TF_hits)    
            
    
    
    ##overlap hits
    if level == 'genes':
        overlap = TF_hits.intersection(GRE_hits)
        k_both = len(overlap)
        hpd = ss.hypergeom(N, k_TF, k_GRE) #hypergeometric distribution of drawing k TFbound genes (total number k_TF) with k_GRE draws given N total genes
        pv = hpd.sf(k_both-1) #p(k >= kboth); is identical to 1-cdf(kboth-1); sf = survival function = 1-cdf
    else:
        if position_res == 'hits':
            nn = n - N*(l+1) #all the possible hit positions for a GRE of length l in N gene promoters 
            k_both = overlap_positions(GRE_hits, TF_hits) 
            hpd = ss.hypergeom(nn, k_TF, k_GRE)
            pv = hpd.sf(k_both-1)
        else:
            k_both = len(TF_hits.intersection(GRE_hits))
            hpd = ss.hypergeom(n, k_TF, k_GRE)
            pv = hpd.sf(k_both-1)
    return [GRE_ID, TF_ID, k_GRE, k_TF, k_both, pv]
        

####-------------run function for every combination of TF-GRE pairs-----------------------------####

for lev in ['genes', 'positions']:
    if lev == 'genes':
        match_list_positions = [nan]*(len(TF_IDs)*len(GRE_IDs))
        i=0
        proba = pb.ProgressBar()
        for GRE in proba(GRE_IDs):
            for TF in TF_IDs:
                print('GRE {} vs TF {}'.format(GRE, TF))
                match_list_positions[i] = position_map_TF_GRE(GRE, TF, 1e-5, target = 'promoters', level = lev)
                i+=1
        pickle.dump(match_list_positions, open('{}match_list_positions_{}.p'.format(out_dir, lev), 'bw'))
    else:
        for resolution in ['hits', 'nucleotides', 'peak_center']:
            match_list_positions = [nan]*(len(TF_IDs)*len(GRE_IDs))
            i=0
            proba = pb.ProgressBar()
            for GRE in proba(GRE_IDs):
                for TF in TF_IDs:
                    print('GRE {} vs TF {}'.format(GRE, TF))
                    match_list_positions[i] = position_map_TF_GRE(GRE, TF, 1e-5, target = 'promoters', level = lev, position_resolution = resolution)
                    i+=1
            pickle.dump(match_list_positions, open('{}match_list_positions_{}_{}.p'.format(out_dir, lev, resolution), 'bw'))