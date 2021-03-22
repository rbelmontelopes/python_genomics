#!/usr/local/bin/python3
#!/Users/rblopes/python_genomics
"""
make the file executable if use the python3 path in your PC. to find the path run 'which python3' on console. To make the file executable run 'chmod a+x filename.py. Triple quotes also indicate comments 
"""

import sys
sys.path.append('/Users/rblopes/python_genomics')

print('***** READING SEQUENCES *****')
print('')

try:
    f=open("dna2.fasta")
except IOError:
    print("File does not exist")

seqs={} #sets an empty dictionary for the sequences

for line in f:
    #discard newline after the end if any (discard on enter)
    line=line.rstrip()
    #distinguish header from sequence
    if line[0] == '>':
        #split header line on whitespaces
        words=line.split() 
        # 0 is the first element of list words, [1:] slice at position 1 leaving the '>' out
        name=words[0][1:]
        #add extracted names to seqs
        seqs[name]=''
    #if sequence data instead of header
    else:
        seqs[name]=seqs[name]+line.upper()
        #append the sequence to the name entry

f.close()
"""
#print the sequences
for name, seq in seqs.items():
    print(name)
    print(seq)

    print("")
    """
    
#print number of sequences
print("There are %d sequences" % len(seqs.keys()))
print('')

print('***** ANALYSING SEQUENCE SIZE *****')
print('')
#creates list with all sequence lenghts and print maximum and minimum
seqs_lenght=[]
seqs_names=[]
seq_sizes=[]
for name, seq in seqs.items(): #iterate over all seqs
    slen=len(seq) #get the lenght of the seqs
    n=name #create a name object for the seqs
    seqs_lenght.append(slen) #append the lenghts to previously created list
    seqs_names.append(n) #append_names to previouly created list
# iterate over names and lengths to create list with both

for i in range(0, len(seqs)):
    seq_sizes.append([seqs_names[i], seqs_lenght[i]])
    #print(seqs_names[i], "; lenght:", seqs_lenght[i], 'base pairs')
#print('')
#create a dict of seqs lenghts
seq_sizes_dict=dict(zip(seqs_names,seqs_lenght))
#get minimum seq lenght
min_seq_lenght=min(seqs_lenght)
#get seq with min size from index position for min lenght in seqs_lenght list
smaller_seq=seqs_names[seqs_lenght.index(min_seq_lenght)]
#get maximum seq lenght
max_seq_lenght=max(seqs_lenght)
#get seq with max size from index position for max lenght in seqs_lenght list
larger_seq=seqs_names[seqs_lenght.index(max_seq_lenght)]
print("The smaller sequence is %s, with %d base pairs" % (smaller_seq, min_seq_lenght))
print("The larger sequence is %s, with %d base pairs" % (larger_seq, max_seq_lenght))
print('')



###################################################################
print('***** SEARCHING FOR ORFs  *****')
print('')
#Find start ORF in FRAME 1
#creates an empty dict to add data from the ORFs per sequence
Repeats_Positions_1={}


#find all repeats and creates a list with repeats and a list with seq names corresponding to repeats
#Creates a loop over the sequences names and values
for name, seq in seqs.items():
    # creates a dict entry for the sequence names (as a dict)
    Repeats_Positions_1[name]={}
    #creates a dict entry for a given sequence for the codon1 data (as a dict)
    Repeats_Positions_1[name]['codon1']={}
    #creates a list entry for a given sequence for the codon1 data start positions
    Repeats_Positions_1[name]['codon1']['start']=[]
    #creates a dict entry for a given sequence for the codon1 data stops (as a dict)
    Repeats_Positions_1[name]['codon1']['stop']={}
    #creates a list entry for a given sequence for the codon1 data TAA stop positions
    Repeats_Positions_1[name]['codon1']['stop']['TAA']=[]
    #creates a list entry for a given sequence for the codon1 data TAG stop positions
    Repeats_Positions_1[name]['codon1']['stop']['TAG']=[]
    #creates a list entry for a given sequence for the codon2 data TAG stop positions
    Repeats_Positions_1[name]['codon1']['stop']['TGA']=[]
    #creates a dict entry for a given sequence for the codon2 data (as a dict) (the rest as for codon1)
    Repeats_Positions_1[name]['codon2']={}
    Repeats_Positions_1[name]['codon2']['start']=[]
    Repeats_Positions_1[name]['codon2']['stop']={}
    Repeats_Positions_1[name]['codon2']['stop']['TAA']=[]
    Repeats_Positions_1[name]['codon2']['stop']['TAG']=[]
    Repeats_Positions_1[name]['codon2']['stop']['TGA']=[]
    #creates a dict entry for a given sequence for the codon3 data (as a dict)
    Repeats_Positions_1[name]['codon3']={}
    Repeats_Positions_1[name]['codon3']['start']=[]
    Repeats_Positions_1[name]['codon3']['stop']={}
    Repeats_Positions_1[name]['codon3']['stop']['TAA']=[]
    Repeats_Positions_1[name]['codon3']['stop']['TAG']=[]
    Repeats_Positions_1[name]['codon3']['stop']['TGA']=[]

    # start in codon 3
    # get first start position
    pos=seqs[name].find('ATG', 0)
    # if first start position % 3 == 2 if belongs to codon3 and will be appended to the start list in the dict
    if pos % 3 == 2:
        p_ = pos +1 #to account for 0 index
        Repeats_Positions_1[name]['codon3']['start'].append(p_)
    # creates a while loop to continue to check for start positions
    while pos >-1:
        pos=seqs[name].find('ATG',pos+3)
        if pos > 0:
            p = pos+1 #to account for 0 index
            pT = pos%3 # filter only for codon 3 positions
            if pT == 2: # filter only for codon 3 positions
                Repeats_Positions_1[name]['codon3']['start'].append(p) # append the position to the start list in the dict
            else:
                pass
    # repeat same steps from above for TAA stops
    pos1=seqs[name].find('TAA', 0)
    if pos1 % 3 == 2:
        p_1 = pos1 +3 # +3 to acount for the last position (pos1 will be the T, and pos1+3 == last A)
        Repeats_Positions_1[name]['codon3']['stop']['TAA'].append(p_1)
    while pos1 >0:
        pos1=seqs[name].find('TAA',pos1+3) 
        if pos1 > 0:
            p1 = pos1+3
            pT1 = pos1%3  
            if pT1 == 2:
                Repeats_Positions_1[name]['codon3']['stop']['TAA'].append(p1)

    pos2=seqs[name].find('TAG', 0)
    if pos2 % 3 == 2:
        p_2 = pos2 +3
        Repeats_Positions_1[name]['codon3']['stop']['TAA'].append(p_2)
    while pos2 >-1:
        pos2=seqs[name].find('TAG',pos2+3)
        if pos2 > 0:
            p2 = pos2+3
            pT2 = (pos2)%3 
            if pT2 == 2:
                Repeats_Positions_1[name]['codon3']['stop']['TAG'].append(p2)

    pos3=seqs[name].find('TGA', 0)
    if pos3 % 3 == 2:
        p_3 = pos3 +3
        Repeats_Positions_1[name]['codon3']['stop']['TGA'].append(p_3)
    while pos3 >-1:
        pos3=seqs[name].find('TGA',pos3+3)
        if pos3 > 0:
            p3 = pos3+3
            pT3 = (pos3)%3
            if pT3 == 2:
                Repeats_Positions_1[name]['codon3']['stop']['TGA'].append(p3)

###Codon 2. Same that for codon 3, but filter is position % 3 == 1

    pos4=seqs[name].find('ATG', 1)
    if pos4 % 3 == 1:
        p_4 = pos4 +1
        Repeats_Positions_1[name]['codon2']['start'].append(p_4)
    while pos4 >-1:
        pos4=seqs[name].find('ATG',pos4+3)
        if pos4 > 0:
            p4 = pos4+1
            pT4 = pos4%3
            if pT4 == 1: # filter only for codon 2 positions
                Repeats_Positions_1[name]['codon2']['start'].append(p4)
    
    pos5=seqs[name].find('TAA', 1)
    if pos5 % 3 == 1:
        p_5 = pos5 +3
        Repeats_Positions_1[name]['codon2']['stop']['TAA'].append(p_5)
    while pos5 >-1:
        pos5=seqs[name].find('TAA',pos5+3)
        if pos5 > 0:
            p5 = pos5+3
            pT5 = p5%3
            if pT5 == 1:
                Repeats_Positions_1[name]['codon2']['stop']['TAA'].append(p5)

    pos6=seqs[name].find('TAG', 1)
    if pos6 % 3 == 1:
        p_6 = pos6 +3
        Repeats_Positions_1[name]['codon2']['stop']['TAG'].append(p_6)
    while pos6 >-1:
        pos6=seqs[name].find('TAG',pos6+3)
        if pos6 > 0:
            p6 = pos6+3
            pT6 = p6%3
            if pT6 == 1:
                Repeats_Positions_1[name]['codon2']['stop']['TAG'].append(p6)

    pos7=seqs[name].find('TGA', 1)
    if pos7 % 3 == 1:
        p_7 = pos7 +3
        Repeats_Positions_1[name]['codon2']['stop']['TGA'].append(p_7)
    while pos7 >-1:
        pos7=seqs[name].find('TGA',pos7+3)
        if pos7 > 0:
            p7 = pos7+3
            pT7 = p7%3
            if pT7 == 1:
                Repeats_Positions_1[name]['codon2']['stop']['TGA'].append(p7)

###Codon 1. Same that for codon 3, but filter is position % 3 == 0

    pos8=seqs[name].find('ATG', 2)
    if pos8 % 3 == 0:
        p_8 = pos8 +1
        Repeats_Positions_1[name]['codon1']['start'].append(p_8)
    while pos8 >-1:
        pos8=seqs[name].find('ATG',pos8+3)
        if pos8 > 0:
            p8 = pos8+1
            pT8 = pos8%3
            if pT8 == 0: # filter only for codon 1 positions
                Repeats_Positions_1[name]['codon1']['start'].append(p8)
    
    pos9=seqs[name].find('TAA', 2)
    if pos9 % 3 == 0:
        p_9 = pos9 +3
        Repeats_Positions_1[name]['codon1']['stop']['TAA'].append(p_9)
    while pos9 >-1:
        pos9=seqs[name].find('TAA',pos9+3)
        if pos9 > 0:
            p9 = pos9+3
            pT9 = p9%3
            if pT9 == 0:
                Repeats_Positions_1[name]['codon1']['stop']['TAA'].append(p9)

    pos10=seqs[name].find('TAG', 2)
    if pos10 % 3 == 0:
        p_10 = pos10 +3
        Repeats_Positions_1[name]['codon1']['stop']['TAG'].append(p_10)
    while pos10 >-1:
        pos10=seqs[name].find('TAG',pos10+3)
        if pos6 > 0:
            p10 = pos10+3
            pT10 = p10%3
            if pT10 == 0:
                Repeats_Positions_1[name]['codon1']['stop']['TAG'].append(p10)

    pos11=seqs[name].find('TGA', 2)
    if pos11 % 3 == 0:
        p_11 = pos11 +3
        Repeats_Positions_1[name]['codon1']['stop']['TGA'].append(p_11)
    while pos11 >-1:
        pos11=seqs[name].find('TGA',pos11+3)
        if pos11 > 0:
            p11 = pos11+3
            pT11 = p11%3
            if pT11 == 0:
                Repeats_Positions_1[name]['codon1']['stop']['TGA'].append(p11)


#calculate OFR lenghts
print('************ ORFs in CODON 1 ************')
# creates a general list to be populated with all ORF sizes for all sequences
ORF_sizes_codon1=[]
#generates a loop to get the values from the dict created with all starts and stops per sequence
for name, codon in Repeats_Positions_1.items():
    #generates a list with all starts for a given sequence (resets at each sequence)
    start1=Repeats_Positions_1[name]['codon1']['start']
    #generates a list with all stops for a given sequence (resets at each sequence)
    stop1=sorted(Repeats_Positions_1[name]['codon1']['stop']['TAA']+Repeats_Positions_1[name]['codon1']['stop']['TAG']+Repeats_Positions_1[name]['codon1']['stop']['TGA'])
    # create a entry in the dict to add the ORF start position, stop position and lenght
    Repeats_Positions_1[name]['codon1']['ORF lenghts']=[]
    # create a list to append all values of ORF start, stop, and lenght for a given sequence
    ORF_len1=[]

    
    #Analyse the ORFs read lenght for codon1
    print("")
    print('Sequence: ', name)
    print("                       ____________codon1____________")
    
    # if no start or stop ORF skip to next sequence
    if start1 == [] or stop1 == []: 
        print('No ORF start or end positions')
        continue
    else:
        pass
    
    # if min start larger than max stop ORF skip to next sequence
    if min(start1) > max(stop1): 
        print('No ORF start or end positions')
        continue
    # else generate a dict with starts and stops for given sequence
    else:
        k=dict.fromkeys(stop1,start1) 
    # initiate a list for start positions which will be reinitiated in the next sequence
    start=[]
    # initiate a list for stop positions which will be reinitiated in the next sequence
    stop=[]
    # get min stop ORF position for adding later
    min_Key=min(k.keys())
    # get max stop ORF position for adding later
    max_Key=max(k.keys())  

    # iterate the dict with starts and stops
    for key, values in k.items():
         # if the stop value is smaller than the first start position, skip, else continue
        if key < values[0]:
            pass
        else:
            # creates loop to iterate over all start values
            for i in range(0, len(values)):
                # if the stop key is smaller than the start value, get index of value and slice the start positions list using the index and append values to start and stop lists
                if key < values[i]: 
                    idx_val_min=values.index(values[i])
                    val_min=values[:idx_val_min]
                    start.append(val_min)
                    stop.append(key)
                    break
# if the min stop key position was not added above, this will add the key if its larger than the max start value            
    if min_Key > k[min_Key][-1]: 
        start.append(k[min_Key])
        stop.append(min_Key)
    else:
        pass
# if the max stop key was not added above, this will add the key if its larger than the max start value    
    if  max_Key > k[max_Key][-1]: 
        start.append(k[max_Key])
        stop.append(max_Key)
    else:
        pass
# generates a new list with first filtering of stops and starts            
    almost=[]
    for i in range(len(stop)):
        alm=stop[i],start[i]
        almost.append(alm)
        almost=sorted(almost)
        
# filter the list from above and generates a new list with only the valid starts and stops
    closer=[]
    # check if first sublist of starts and stops fits the criteria (stop being larger than start) cause it is not contemplated in the loop below, and if so add to the new list
    if almost[0][0] > almost[0][1][-1]:
                clr=almost[0][0] , almost[0][1]
                closer.append(clr)
    #creates a loop to check for values of start that are larger than the the anterior stop, get the index for the value, slice the start values to include only values larger than the anterior stop position,  append start values and stop to the new list and break to not add the values several times 
    for i in range(0, len(almost)):
        for j in range(len(almost[i][1])):
            if almost[i-1][0] < almost[i][1][j]:
                idx_alm=almost[i][1].index(almost[i][1][j])
                right_slice=almost[i][1][idx_alm:]
                clr= almost[i][0], right_slice
                closer.append(clr)
                break

# creates a loop to calculate the ORF values from the above list                    
    for i in range(0, len(closer)):
            O_len1=[]
            O_len1=(closer[i][0] - min(closer[i][1]))+1 # add one to compensate for zero indexing
            O_len1_start=[]
            O_len1_start=min(closer[i][1]) #get the start position
            ORF_len1_stop=[]
            ORF_len1_stop=closer[i][0] #get the stop position
            O_len1_final=[]
            O_len1_final=O_len1_start, ORF_len1_stop, O_len1 # creates a list with start, stop and ORF lenght
            ORF_len1=[]
            ORF_len1.append(O_len1_final) # appends the list with start, stop and ORF lenght to a list
            Repeats_Positions_1[name]['codon1']['ORF lenghts'].append(ORF_len1) #add start, stop and ORF lenght to the main dict 
            print('ORF start: ', O_len1_start, ' ORF end: ', closer[i][0], ' lenght: ',O_len1, ' base pairs')  #print the values for the sequence
            ORF_s_c1=name, O_len1_start, O_len1 # creates a sublist including sequence name, start position and ORF lenght
            ORF_sizes_codon1.append(ORF_s_c1) #append the sublist above to a general list (which will include all sequences) to be used to calculate the maximum ORF lenght for codon1

#generates a list to include all ORF lenghts for codon1
sizes=[]
#create a loop to append only the ORF lenghts to the list above
for i in range(0, len(ORF_sizes_codon1)):
    sizes.append(ORF_sizes_codon1[i][2])
#get the maximum value from the populated ORF lenghts list 
MAX_ORF1=max(sizes)
#get the index of the maximum value from above
MAX_ORF1_idx=sizes.index(MAX_ORF1)
#print the slice the general list of ORF sizes for codon1 with the maximum value index (as lists are immutable, the index val from sizes list will correspond to the one of the general list)
print(ORF_sizes_codon1[MAX_ORF1_idx])
print('')

                  
#########CODON2
print('************ ORFs in CODON 2 ************')
ORF_sizes_codon2=[]
for name, codon in Repeats_Positions_1.items():

    
    
    start2=Repeats_Positions_1[name]['codon2']['start']    
    stop2=sorted(Repeats_Positions_1[name]['codon2']['stop']['TAA']+Repeats_Positions_1[name]['codon2']['stop']['TAG']+Repeats_Positions_1[name]['codon2']['stop']['TGA'])
    Repeats_Positions_1[name]['codon2']['ORF lenghts']=[]

    ORF_len2=[]

    print('')
    print('Sequence: ', name)
    print("                       ____________codon2____________")

    if start2 == [] or stop2 == []:
        print('No ORF start or end positions')
        continue
    else:
        pass
    
    if min(start2) > max(stop2):
        print('No ORF start or end positions')
        continue
    else:
        k=dict.fromkeys(stop2,start2)

    start=[]
    stop=[]
    min_Key=min(k.keys())
    max_Key=max(k.keys())

    for key, values in k.items():
        if key < values[0]:
            continue
        else: 
            for i in range(0, len(values)):
                if key < values[i]:
                    idx_val_min=values.index(values[i])
                    val_min=values[:idx_val_min]
                    start.append(val_min)
                    stop.append(key)
                    break
            
    if min_Key > k[min_Key][-1]:
        start.append(k[min_Key])
        stop.append(min_Key)
    else:
        pass
    
    if  max_Key > k[max_Key][-1]:
        start.append(k[max_Key])
        stop.append(max_Key)
    else:
        pass
                
    almost=[]
    for i in range(len(stop)):
        alm=stop[i],start[i]
        almost.append(alm)
        almost=sorted(almost)
        

    closer=[]
    if almost[0][0] > almost[0][1][-1]:
                clr=almost[0][0] , almost[0][1]
                closer.append(clr)
    for i in range(0, len(almost)):
        for j in range(len(almost[i][1])):
            if almost[i-1][0] < almost[i][1][j]:
                idx_alm=almost[i][1].index(almost[i][1][j])
                right_slice=almost[i][1][idx_alm:]
                clr= almost[i][0], right_slice
                closer.append(clr)
                break

                
       
            
                
    for i in range(0, len(closer)):
        if closer[i][1] == []:
            pass
        else:
            O_len2=[]
            O_len2=(closer[i][0] - min(closer[i][1]))+1
            O_len2_start=[]
            O_len2_start=min(closer[i][1])
            ORF_len2_stop=[]
            ORF_len2_stop=closer[i][0]
            O_len2_final=[]
            O_len2_final=O_len2_start, ORF_len2_stop, O_len2
            ORF_len2=[]
            ORF_len2.append(O_len2_final)
            Repeats_Positions_1[name]['codon2']['ORF lenghts'].append(ORF_len2)
            print('ORF start: ', O_len2_start, ' ORF end: ', closer[i][0], ' lenght: ',O_len2, ' base pairs')
            ORF_s_c2=name, O_len2_start, O_len2
            ORF_sizes_codon2.append(ORF_s_c2)
sizes=[]
for i in range(0, len(ORF_sizes_codon2)):
    sizes.append(ORF_sizes_codon2[i][2])
MAX_ORF2=max(sizes)
MAX_ORF2_idx=sizes.index(MAX_ORF2)
print(ORF_sizes_codon2[MAX_ORF2_idx])
print('')

#########CODON3
print('')
print('************ ORFs in CODON 3 ************')
ORF_sizes_codon3=[]
for name, codon in Repeats_Positions_1.items():
    start3=Repeats_Positions_1[name]['codon3']['start']    
    stop3=sorted(Repeats_Positions_1[name]['codon3']['stop']['TAA']+Repeats_Positions_1[name]['codon3']['stop']['TAG']+Repeats_Positions_1[name]['codon3']['stop']['TGA'])
    Repeats_Positions_1[name]['codon3']['ORF lenghts']=[]
    Repeats_Positions_1[name]['codon3']['Max ORF lenght']=[]

    ORF_len3=[]
    print('')
    print('Sequence: ', name)
    print('')
    print("                       ____________codon3____________")

    
    if start3 == [] or stop3 == []:
        print('No ORF start or end positions')
        continue
    else:
        pass
    
    if min(start3) > max(stop3):
        print('No ORF start or end positions')
        continue
    else:
        k=dict.fromkeys(stop3,start3)

    start=[]
    stop=[]
    min_Key=min(k.keys())
    max_Key=max(k.keys())

    for key, values in k.items():
        if key < values[0]:
            pass
        else: 
            for i in range(0, len(values)):
                if key < values[i]:
                    idx_val_min=values.index(values[i])
                    val_min=values[:idx_val_min]
                    start.append(val_min)
                    stop.append(key)
                    break
            
    if min_Key > k[min_Key][-1]:
        start.append(k[min_Key])
        stop.append(min_Key)
    else:
        pass
    
    if  max_Key > k[max_Key][-1]:
        start.append(k[max_Key])
        stop.append(max_Key)
    else:
        pass
                
    almost=[]
    for i in range(len(stop)):
        alm=stop[i],start[i]
        almost.append(alm)
        almost=sorted(almost)
        

    closer=[]
    if almost[0][0] > almost[0][1][-1]:
                clr=almost[0][0] , almost[0][1]
                closer.append(clr)
    for i in range(0, len(almost)):
        for j in range(len(almost[i][1])):
            if almost[i-1][0] < almost[i][1][j]:
                idx_alm=almost[i][1].index(almost[i][1][j])
                right_slice=almost[i][1][idx_alm:]
                clr= almost[i][0], right_slice
                closer.append(clr)
                break
                
    for i in range(0, len(closer)):
        if closer[i][1] == []:
            pass
        else:
            O_len3=[]
            O_len3=(closer[i][0] - min(closer[i][1]))+1
            O_len3_start=[]
            O_len3_start=min(closer[i][1])
            ORF_len3_stop=[]
            ORF_len3_stop=closer[i][0]
            O_len3_final=[]
            O_len3_final=O_len3_start, ORF_len3_stop, O_len3
            ORF_len3=[]
            ORF_len3.append(O_len3_final)
            Repeats_Positions_1[name]['codon3']['ORF lenghts'].append(ORF_len3)
            print('ORF start: ', O_len3_start, ' ORF end: ', closer[i][0], ' lenght: ',O_len3, ' base pairs')
            ORF_s_c3=name, O_len3_start, O_len3 
            ORF_sizes_codon3.append(ORF_s_c3) #creates list with all orfs for codon3, with name and start position

#code to find the largest of in the list of ORF in codon3 created above
sizes=[]
for i in range(0, len(ORF_sizes_codon3)):
    sizes.append(ORF_sizes_codon3[i][2])
MAX_ORF3=max(sizes)
MAX_ORF3_idx=sizes.index(MAX_ORF3)
print(ORF_sizes_codon3[MAX_ORF3_idx])
print('')
            
    
############################################################################          

#lenght of the tandem to search
repeat_num=7 #6, 12, 7
print('')
print('***** SEARCHING FOR REPEATS WITH %d NUCLEOTIDES *****' % repeat_num)

#Count repeats of n NUCLEOTIDES for the file as a whole (all sequences)
#creates empty tandem dict
tandem_dic = {}
#creates empty list for tandems
tandem_list=[]
#creates empty list for tandem sequence names
tandem_names=[]


#creates a loop to iterate over all sequences
for name, seq in seqs.items():
    #creates an entry for each sequence as a subdict in the tandem dict
    tandem_dic[name]={}
    #creates a list to append the tandems in each sequence
    tandem_dic[name]['tandem']=[]

    #loops over all seqs and add tandems to a general list (not associated with sequence or position)
    for i in range((len(seq)-repeat_num)):  
        tandem=(seq[i:(i+repeat_num)])
        tandem_list.append(tandem)

    #loop over a sequence and add tandems to the tandem dict. Still need to be correct
    tandem_names=[] # to reset tandem names list for each seq
    for i in range((len(seqs[name])-repeat_num)):  
        tandem=(seq[i:(i+repeat_num)])
        tandem_names.append(tandem)
        tandem_dic[name]['tandem'].append(tandem_names)
        
#creates a dict of repeats with repeats as keys and count as values
repeat_counter={}
for repeat in tandem_list: 
    if repeat in repeat_counter:
        repeat_counter[repeat] +=1
    else:
        repeat_counter[repeat] =1
#sort most common firt            
most_common = sorted(repeat_counter, key = repeat_counter.get, reverse=True)
#slice the 6 most common repeats
top_6 = most_common[0:6]
#print repeat name + count
for item in (top_6):
    print(item, repeat_counter[item])       







"""    
  
general_dict={}
for name, seq in seqs.items():
    general_dict[name]={}
    general_dict[name].update(sequence=seq)
    general_dict[name].update(seq_lenght=seq_sizes_dict[name]) 
    general_dict[name].update(start_ORFS=start_ORFs[name])
    general_dict[name].update(TAA_stops=end_ORFstTAA[name])
    general_dict[name].update(TAG_stops=end_ORFstTAG[name])
    general_dict[name].update(TGA_stops=end_ORFstTGA[name])


"""
