import os, subprocess
from progressbar import ProgressBar
pbar = ProgressBar()
import pickle


tomtom_path = '/Users/rkoch/Programs/meme/bin/tomtom'
#Here the GRE MEME output files are stored
GRE_directory = '/Users/rkoch/Documents/Data_and_Scripts/motif_clusters_24/'
#Here the TF MEME output files are stored
TF_directory = '/Users/rkoch/Documents/Data_and_Scripts/MEME_OUTPUT/'
#this is where the output will be stored
out_dir = '/Users/rkoch/Documents/Data_and_Scripts/out/match_lists/'

nan = float('NaN')

#helper function
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
###---------------------------------main-function I----------------------------------------------###     
def run_tomtom(TF_file, GRE_file, e_thresh, dist_method,
               min_overlap, TF_ID= 'all', extended_info = True):
    
    """
    a wrapper around the tomtom script
    
    input: 
        TF_file - a file containing two TF chip motifs (MEME output)
        GRE_file - a file containing a GRE motif (MEME output)
        e_thresh - a threshold for the E-value of motif match
        dist_method - the method used to calculate the distance between the sequences 
        min_overlap - the minimum overlap between the sequences
        TF_ID - either 1, 2 or 'all', depending on which motif(s) is/are used. 1 and 2 are used when run_tomtom is called from tomtom_splicer
        extended info - boolean, should meta-info w.r.t. motifs be added to the output?
        
    output: 
        a list of strings with each string corresponding to the output of a tomtom comparison TF_query vs GRE_query : 
            TF_motif_ID, 
            GRE_motif_ID(always1), 
            offset, 
            pval, 
            qval, 
            eval, 
            overlap, 
            TF sequence, 
            GRE sequence, 
            orientation,
            if extendex_info == True: 
                Evalue TF motif, 
                number of seqs used to create TF motif, 
                Evalue GRE motif, 
                number of ses used to create GRE motif
            
    """


    if TF_ID not in [1,2,'all']:
        print("ERROR: TF_ID must be 1, 2, or 'all'")
        return
     
    
    if extended_info:
        TF_content = []
        with open(TF_file) as f:
            TF_content = f.readlines()
        l_o_i = [x.split() for x in TF_content if x.startswith('MOTIF')]
        TF_evals = [float(x[-1]) for x in l_o_i]
        TF_number_of_seqs = [int(x[7]) for x in l_o_i]
        with open(GRE_file) as f:
            GRE_content = f.readlines()
        l_o_i = [x.split() for x in GRE_content if x.startswith('MOTIF')]
        GRE_eval = float(l_o_i[0][-1])
        GRE_number_of_seqs = int(l_o_i[0][8])
    
    

    #two variants, one scanning only one motif, the other one all motifs
    if TF_ID in [1,2]:
        command = [tomtom_path,
                   '-no-ssc', 
                   '-m', str(TF_ID),
                   '-oc', '.',
                   '-verbosity', '1',
                   '-evalue',
                   '-thresh', '%.3f' % e_thresh,        
                   '-dist', dist_method,
                   '-min-overlap', '%d' % min_overlap,
                   '-text',
                   TF_file,GRE_file]
    else:
        command = [tomtom_path,
                   '-no-ssc', 
                   '-oc', '.',
                   '-verbosity', '1',
                   '-evalue',
                   '-thresh', '%.3f' % e_thresh,        
                   '-dist', dist_method,
                   '-min-overlap', '%d' % min_overlap,
                   '-text',
                   TF_file,GRE_file]
            

    #collect tomtom output
    output = subprocess.check_output(command).decode('utf-8')

    #reformat output before return
    lines = output.split('\n')[1:-1]
    if not lines:
        return None
    elif not extended_info:
        return lines 
    else:
        for l, eva, no_seq in zip(range(len(lines)), TF_evals, TF_number_of_seqs):
            lines[l] = lines[l] + '\t{}\t{}\t{}\t{}'.format(eva, no_seq, GRE_eval, GRE_number_of_seqs)
        return lines


###---------------------------------main-function I----------------------------------------------###     
def tomtom_splicer(TF_query,GRE_query,variant = "both"): 
    
    """ 
    this function runs tomtom on a GRE and a TF motif, truncates the sequences based on the naive
    run and reruns tomtom with the truncated sequence
    
    input:
        TF_query - path of a TF MEME-file-containing folder
        GRE_query - path of a GRE MEME-file
        variant - one of "GRE" or "both", indicating whether both motifs or only the GRE should be cut
        
    output:
        a list of strings with each string corresponding to the output of a tomtom comparison TF_query vs spliced GRE_query
    
    """  
    output = []


    #produce 'regular' tomtom output
    hi = run_tomtom(TF_query, GRE_query,  1000, 'pearson', 5, TF_ID = 'all', extended_info = False)
    
    if hi == None:
        return None
    #iterate over both result rows, i.e. both TF motifs
    for ind in range(len(hi)):
        row = hi[ind].split()
        if len(row) != 10:  #test for acceptable format
            return None
        ident = int(row[0]) #query ID = TF motif ID
        off = int(row[2])   #offset
        o = int(row[6])     #overlap
        m = len(row[7])     #length of query consensus = TF motif consensus
        n = len(row[8])     #length target consensus = GRE consensus
        sign = row[9]       #orientation: plus or minus strand
        if sign == '-':     
            revcomp = True
        else:
            revcomp = False
        
        #load the GRE input file for subsequent manipulation based on native tomtom output
        with open(GRE_query, newline='') as inputfile2:
            GRE_data = inputfile2.readlines()
        #anchor line indices for PSSM data
        indis_GRE = [i for i, j in enumerate(GRE_data) if ('position-specific' in j and 'matrix' in j)]
        #to start from bottom of file so that integrity is conserved
        indis_GRE.reverse()
        
        first_G = True
        
        
        if variant == "both":
            #load the TF input file for subsequent manipulation based on native tomtom output
            with open(TF_query) as inputfile1:
                TF_data = inputfile1.readlines()
            indis_TF = [i for i, j in enumerate(TF_data) if ('position-specific' in j and 'matrix' in j and str(ident) in j)]
            indis_TF.reverse()
            first_T = True
        
         
        
        ############################ Case 1 ############################

        #           TF_______________ m
        #               GRE________   n
        if o == n:
            print('Case I')
            case = 1
            if variant == 'both':
                for k in indis_TF:    
                    #substitute descriptive lines in file
                    TF_data.pop(k + 2)
                    if first_T:
                        TF_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                        first_T = False
                    else:   
                        TF_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    #remove non-overlapping positions
                    for i in range(m-n-abs(off)): # need to remove #n-m-off elements on the right side
                        TF_data.pop(k + 3 + abs(off) + n)
                    for i in range(abs(off)): # need to remove #off elements on the left side
                        TF_data.pop(k + 3)
                

        ############################ Case 2 ############################
        
        #          ....off..TF___________.n-off-m..        m
        #       GRE________________________________        n
        
        if o == m: 
            print('Case II')
            case = 2
            for k in indis_GRE:    
                #substitute descriptive lines in file
                GRE_data.pop(k + 2)
                if first_G:
                    GRE_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    first_G = False
                else:   
                    GRE_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                #remove non-overlapping positions
                if not revcomp:
                    for i in range(n-m-off): # need to remove #n-m-off elements on the right side
                        GRE_data.pop(k + 3 + off + m)
                    for i in range(off): # need to remove #off elements on the left side
                        GRE_data.pop(k + 3)
                if revcomp:
                    for i in range(off):
                        GRE_data.pop(k + 3 + n - off) # need to remove #off elements on the right side
                    for i in range(n-m-off):
                        GRE_data.pop(k + 3) # need to remove #n-m-off elements on the left side
                        
                        
                        
                        
        ############################ Case 3 ############################
        
        #           TF_______________                          m
        #               GRE_____________________________       n
            
        if off < 0 and o < n: 
            print('Case III')
            case = 3
            for k in indis_GRE:                
                #substitute descriptive lines in file
                GRE_data.pop(k + 2)
                if first_G:
                    GRE_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    first_G = False
                else:   
                    GRE_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                #remove non-overlapping positions
                if not revcomp:
                    for i in range(n-o): #all the positions that need to be removed from the end of GRE
                        GRE_data.pop(k + 3 + o)
                if revcomp:
                    for i in range(n-o):
                        GRE_data.pop(k + 3)
            if variant == 'both':
                for k in indis_TF:
                    TF_data.pop(k+2)
                    if first_T:
                        TF_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                        first_T = False
                    else:   
                        TF_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    for i in range(abs(off)):
                        TF_data.pop(k + 3)
                    
           
        
                
 
   
        ############################ Case 4 ############################
    
        #                        off            TF_______________   m
        #               GRE______________________________           n
        
        if off >= 0 and off + m > n:
            print('Case IV')
            case = 4
            for k in indis_GRE:                
                #substitute descriptive lines in file
                GRE_data.pop(k + 2)
                if first_G:
                    GRE_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    first_G = False
                else:   
                    GRE_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                #remove non-overlapping positions
                if not revcomp:
                    for i in range(off):
                        GRE_data.pop(k + 3)
                if revcomp:
                    for i in range(off):
                        GRE_data.pop(k + 3 + o)
            if variant == 'both':
                for k in indis_TF:
                    #substitute descriptive lines in file
                    TF_data.pop(k + 2)
                    if first_T:
                        TF_data.insert(k + 2, 'letter-probability matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))                         
                        first_T = False
                    else:
                        TF_data.insert(k + 2, 'log-odds matrix: alength= 4 w= %s n= 1658 bayes= 5.8633 E= 2.1e-118 \n' %(o))
                    #remove non-overlapping positions
                    for i in range(m-o):
                        TF_data.pop(k + 3 + o)
                
                                    
        with open("{}truncated_GRE_data_{}.txt".format(out_dir,ident), "w") as f:
           for row in GRE_data:
                f.write(row)
        
        with open("{}truncated_TF_data_{}.txt".format(out_dir,ident), 'w') as f:
            for row in TF_data:
                f.write(row)
        
        #print(ident)    
        runling = run_tomtom("{}truncated_TF_data_{}.txt".format(out_dir,ident), '{}truncated_GRE_data_{}.txt'.format(out_dir,ident), 1000, 'pearson', 5, TF_ID = ident, extended_info = False)

        if runling != None:
            runling = runling[0] + '\t{}'.format(case)
            output.append(runling)
              
    return output



####--------run function for all peak sets for every combination of TF-GRE pairs----------------####


#the sets of peaks on which motif generation is based
types =  ['out', 'in_and_expressed', 'all','in']
#load the GRE motifs 
targets = os.listdir(GRE_directory) 
targets = [x for x in targets if x.endswith('memeOut.txt')]
#load the palindromic&nonpalindromic TF motifs from TFOE data
for typ in types:
    print(typ)
    #TF motifs specific for the peak set 'typ'
    queries = os.listdir('{}{}'.format(TF_directory,typ)) 
    queries = [x for x in queries if "Rv" in x]
    
    
    #create list of matches without cutting any motifs
    pbar = ProgressBar()
    match_list = [nan]* (len(targets)*len(queries)*2)
    i = 0
    for GRE in pbar(targets):
        for TF in queries:
            print('{} vs {} naive tomtom'.format(GRE, TF))
            GRE_ID = int(GRE.split("_")[0])
            TF_ID, shape = TF.split('_')
            head = [GRE_ID,TF_ID, shape]
            match = run_tomtom('{}{}/{}/meme.txt'.format(TF_directory, typ, TF), '{}{}'.format(GRE_directory,GRE), 1000, 'pearson', 5, extended_info = True)
            #transfer tomtom output into a list of outputs
            if match is not None:
                for line in match:
                    line = line.split('\t')
                    for x in range(len(line)):       #next ~12 lines: cleaning & variable formatting
                        y = line[x]
                        if is_number(y):
                            if y.isnumeric():
                                line[x] = int(y)
                            else:
                                line[x] = float(y)
                                if line[x] < 0:
                                    line[x] = int(line[x])
                        elif y not in ('-', '+'):
                            line[x] = len(y)
                    line.pop(1)                                                 #GRE ID, 1 in any case -> discard
                    if line[8] == '+':                  
                        center = line[1] + line[6]/2 + 0.5                     #center of binding wrt GRE: off + len(TF)/2  + 0.5
                    if line[8] == '-':
                        center = line[7] + 1 - (line[1] + line[6]/2 + 0.5)    #center of binding: len(GRE)+1 - (off + len(TF)/2  + 0.5)
                    insert = head + line + [center]
                    match_list[i] = insert
                    i += 1
            else:
                for j in [1,2]:
                    match_list[i] = head+[j]+[nan]*13
                    i += 1
    
    
    #to store the match list
    #match_list = [x for x in match_list if type(x) == list]            
    pickle.dump(match_list, open('{}/full_match_list_{}.p'.format(out_dir,typ), 'bw'))   
    #to restore the match list
    #match_list = pickle.load(open('{}/full_match_list_{}.p'.format(out_dir,typ), 'rb'))

   
    #create list of matches after cutting both motifs
    vari = 'both'
    pbar = ProgressBar()            
    spliced_match_list = [nan]*(len(targets)*len(queries)*2)
    i = 0
    for GRE in pbar(targets):
        for TF in queries:
            print('{} vs {} spliced tomtom, variant = {}'.format(GRE, TF, vari))
            GRE_ID = int(GRE.split("_")[0])
            TF_ID, shape = TF.split('_')
            head = [GRE_ID,TF_ID, shape]
            match = tomtom_splicer(TF_query = '{}{}/{}/meme.txt'.format(TF_directory,typ,TF), GRE_query = '{}{}'.format(GRE_directory, GRE), variant = vari) 
            if match is not None:
                for line in match:
                    line = line.split('\t')
                    for x in range(len(line)):       #next ~12 lines: cleaning & variable formatting
                        y = line[x]
                        if is_number(y):
                            if y.isnumeric():
                                line[x] = int(y)
                            else:
                                line[x] = float(y)  
                                if line[x] < 0:
                                    line[x] = int(line[x])
                        elif y not in ('-', '+'):
                            line[x] = len(y)
                    line.pop(1)                     
                    insert = head + line
                    spliced_match_list[i] = insert
                    i += 1
            else:
                for j in [1,2]:
                    spliced_match_list[i] = head+[j]+[nan]*9
                    i += 1
    spliced_match_list_both = [x for x in spliced_match_list if type(x) == list]
    pickle.dump(spliced_match_list, open('{}/full_double_spliced_match_list_{}.p'.format(out_dir,typ), 'bw'))
    
    
    #create list of matches after cutting only GRE motifs
    vari = 'GRE'
    pbar = ProgressBar()            
    spliced_match_list = [nan]*(len(targets)*len(queries)*2)
    i = 0
    for GRE in pbar(targets):
        for TF in queries:
            print('{} vs {} spliced tomtom, variant = {}'.format(GRE, TF, vari))
            GRE_ID = int(GRE.split("_")[0])
            TF_ID, shape = TF.split('_')
            head = [GRE_ID,TF_ID, shape]
            match = tomtom_splicer('{}{}/{}/meme.txt'.format(TF_directory,typ,TF), '{}{}'.format(GRE_directory, GRE), 1000, 5, variant = vari) 
            if match is not None:
                for line in match:
                    line = line.split('\t')
                    for x in range(len(line)):       #next ~12 lines: cleaning & variable formatting
                        y = line[x]
                        if is_number(y):
                            if y.isnumeric():
                                line[x] = int(y)
                            else:
                                line[x] = float(y)  
                                if line[x] < 0:
                                    line[x] = int(line[x])
                        elif y not in ('-', '+'):
                            line[x] = len(y)
                    line.pop(1)                     
                    insert = head + line
                    spliced_match_list[i] = insert
                    i += 1
            else:
                for j in [1,2]:
                    spliced_match_list[i] = head+[j]+[nan]*9
                    i += 1
    spliced_match_list_GRE = [x for x in spliced_match_list if type(x) == list]
    pickle.dump(spliced_match_list, open('{}/full_GRE_spliced_match_list_{}.p'.format(out_dir,typ), 'bw'))
    



    #combining results of nonspliced and spliced matching
    combined_match_list = [nan]*(len(match_list))
    pbar = ProgressBar()
    ix = 0
    print('combining the naive and spliced matches for type {}'.format(typ))
    for i in pbar(match_list):
        for j, k in zip(spliced_match_list_both,spliced_match_list_GRE):
            if (i[0:4] == j[0:4]):
                #entries: GRE_ID, TF ID, TF shape, TF motif ID(i0:4), 
                #case of overlap(j-1), 
                #length overlap(i8),  
                #length of GRE unspliced(i10), 
                #length of GRE spliced(j10), 
                #number of GRE nt spliced(i10-j10), 
                #length of TFmotif unspliced(i9), 
                #length of TFmotif spliced(j9), 
                #number of TF nt spliced(i9-j9),
                #pvalue unspliced(i5), 
                #pvalue spliced both(j5),
                #pvalue spliced GRE only (k5),
                #orientation(i11), 
                #center of binding wrt GRE(i-1), 
                #evalueTF, no_seqsTF, evalueGRE, no_seqsGRE (i12:15)
                combined_match_list[ix] = i[0:4] + [j[-1]] + [i[8]] + [i[10]] + [j[10]] + [i[10]-j[10]]+ [i[9]] + [j[9]] + [i[9]-j[9]] + [i[5]] + [j[5]] +[k[5]] + [i[11]] + [i[-1]] + i[12:16]
                ix += 1
                continue
    combined_match_list = [x for x in combined_match_list if type(x) == list]
    pickle.dump(combined_match_list, open('{}/combined_match_list_{}.p'.format(out_dir, typ), 'bw'))





