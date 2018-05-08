        #!/usr/bin/env python3
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import pickle
#import  csv
#from progressbar import ProgressBar
#pbar = ProgressBar()

Gene = 'Rv0397' #one specific gene 
print(Gene)
INPUTFILE = 'all_genes_fimos/%s.csv' %(Gene) # a table of all nucleotides within Gene, with the number of GREs overlapping significantly at this position
CHIPSEQ_FILE = 'Data/chip_seq_data.csv'


chip_dic = pickle.load(open("Dictionaries/chip_dic.p", "rb" )) ## contains for every TF all positions where chipseq peak was found (+- 12bp)
chip_dict = pickle.load(open("Dictionaries/chip_dict.p", "rb" )) ## contains for every TF all genes where chipseq peaks lies within promoter regulating gene
promoter_dic = pickle.load(open( "Dictionaries/gene_dic.p", "rb" )) ## contains the promoter region for every gene
#pub_chip_dic = pickle.load(open( "Dictionaries/pub_chip_dic.p", "rb" )) ## to load dictionary
#new_ffl_list = pickle.load(open( "Dictionaries/new_ffl_list.p", "rb" )) ## to load dictionary
#data1 = open('Data/de.corems.csv')
#read_data1 = csv.reader(data1)
#data1.close()

def main():
    plot_gres(INPUTFILE, take_top=20, min_count=20)

    



def plot_gres(path, take_top=None, min_count=None):
    # load the GRE-promoterregion overlap data and extract it
    df = pd.read_csv(path, index_col=0)
    gene = sorted(set(df.loc[:,'GENE']))[0] # should be the same along the whole column anyway
    gre_ids = sorted(set(df.loc[:,'id']))
    gre_max = []

    for gre_id in gre_ids:
        pos_and_count = df[df.loc[:,'id'] == gre_id].loc[:,['start','value']]
        y = sorted(list(pos_and_count.loc[:, 'value'].values))
        gre_max.append((max(y), gre_id))

    gre_max = sorted(gre_max, reverse=True)
    if take_top is not None:
        gre_max = gre_max[:take_top]
        top_gre_list = []
        for n in range(len(gre_max)):
            if gre_max[n][0] > min_count:
                top_gre_list.append(gre_max[n][1])
    

    if min_count is not None:
        gre_max = [(y, gre_id) for y, gre_id in gre_max if y >= min_count]
        count_max = [y for y, gre_id in gre_max if y >= min_count]
        count_gres = [gre_id for y, gre_id in gre_max if y >= min_count]
        
#    for gre in count_gres:
#        gre_list.append(gre)
    
    fig = plt.figure()
    ax = plt.subplot(111)

    tss = sorted(set(df.loc[:,'TSS']))[0]
    plot_gre_data(ax, df, map(lambda e: e[1], gre_max), tss)
    rangex = plt.xlim()
    xwindow = (int(tss + rangex[0]), int(tss + rangex[1]))
#    print([tss - 150,tss + 70])

    if Gene[-1] != 'c':
#        print('hi')
        window_start = tss - 150
#        print(window_start)
        window_end = tss + 70
        for chip in chip_dic:
            for pos in chip_dic[chip]:
#                print(pos)
                if (pos[0] + 12) >= window_start and (pos[0] + 12) <= window_end:
                    plt.axvline((pos[0] + 12) - tss)
                    plt.text((pos[0] + 12) - tss,30,chip,rotation=90)
#                    sig_tfs[Gene].append(chip)
    #                    print([chip,pos[0]+12])
                        
    if Gene[-1] == 'c':
#        print('hi')
        window_start = tss - 70
        window_end = tss + 150
        for chip in chip_dic:
            for pos in chip_dic[chip]:
                if (pos[0] + 12) >= window_start and (pos[0] + 12) <= window_end:
                    plt.axvline(tss - (pos[0] + 12))
                    plt.text(tss - (pos[0] + 12),5,chip,rotation=90)
#                    sig_tfs[Gene].append(chip)
    #                    print([chip,pos[0]+12])

#    tf_centers = pd.read_csv(CHIPSEQ_FILE, index_col=0, header=None, names=['center'])
#    tf_centers[(tf_centers['center'] >= xwindow[0]) & (tf_centers['center'] < xwindow[1])]
#    print(xwindow)
#    sig_gres[Gene] = top_gre_list
    #plot_tf_centers(ax, tf_centers)
#    plt.savefig("Fimo Plots/fuzzy_promoters/%s.svg" %(Gene))
    if count_max:
        if max(count_max) > min_count:
            plt.show()
          
def plot_gre_data(ax, df, gre_ids, tss):
    for gre_id in gre_ids:
        pos_and_count = df[df.loc[:,'id'] == gre_id].loc[:,['start','value']]
        x = list(pos_and_count.loc[:, 'start'].values)
        x = list(map(lambda v: v - tss, x))
        y = list(pos_and_count.loc[:, 'value'].values)
        x.append(max(x) + 1)
        y.append(0)
#        print(x,y)
#        ax.plot(x, y, label=gre_id)
        if Gene[-1] != 'c':
            ax.plot(x, y, label=gre_id)
        if Gene[-1] == 'c':
            x = [-c for c in x]
#            y = [-d for d in y]
            ax.plot(x, y, label=gre_id)
#        GRE_counter[gre_id] = GRE_counter[gre_id] + 1

    # position and size the main plot to match up with the legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])

    lgd = ax.legend(ncol=5, loc='upper center', bbox_to_anchor=(0.5, -0.10),
                    fancybox=True, shadow=True, fontsize='x-small')



                #somehow deprecated stuff
# find all positions inside a promoter (target) where TF is supposed to bind based on Chip-Seq data
#def tf_in_promoter(TF,target): 
#    pos_list = []
#    for pos in chip_dic[TF]:
#        if promoter_dic[target][0] <= pos[0] + 12 <= promoter_dic[target][1]:
#            pos_list.append(pos[0] + 12)
#    return pos_list
    
#def plot_tf_centers(ax, tf_centers):
#    ax.text(-140, 210, "Rv0001")
#    ax.plot([-140, -130], [200, 200], label='blabla')

  
            
if __name__ == '__main__':
    main()