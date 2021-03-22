import collections
import matplotlib.pyplot as plt
import random
import statistics
import numpy as np
from scipy import stats
import math
import bisect
import datetime as d
from itertools import permutations


###These functions were created during the Coursera Genomic Data Science Course "Algorithms for DNA sequencing", and some of them are authored by Ben Langmead. Several others were created as part of the course by myself. I tried to indicates the authors to all functions but possibly there are some mismatches in the authors


###func to get the longest common prefix between two seqs
def longestCommonPrefix(s1, s2):
    i=0 # index for the seqs
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
    #while i is smaller than the lenght of the seqs and the characters in s1 and s2 are equal
        i += 1
    return s1[:i] # will return the longest common prefix

###func to check if sequences match
def match (s1, s2):
    if not len(s1) == len(s2):
        return False
    for i in range(len(s1)):
        if not s1[i] == s2[i]:
            return False
    #if the loop above finished without returning false means that all charactes are equal, so return True        
    return True

###func to reverse complement a seq
def reverseComplement(seq):
    complement={'A': 'T', 'T':'A', 'C' : 'G', 'G':'C', 'N':'N'}
    t = ''
    for base in seq:
        t = complement[base] + t # add complement to the begining cause the reversion of the string
    return t

###func to read a genome / fasta file
def readGenome(filename):
    genome=''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>': # if is other line than the first (which starts with the >)
                genome += line.rstrip() #.rstrip trim the end of the line (enter, tab, whitespace)
    return genome

### func to convert chr to ASCII      
def QtoPhred(Q):
    "turn q into Phred+33 ASCII encoded quality"
    return chr(Q + 33) # chr convert integer to ASCII

###func to convert ASCII to char    
def phred33ToQ(qual):
    "Turn phred+33 ASCII-encoded quality into Q"
    return ord(qual)-33 #ORD convert ASCII to integer

###func to read FastQ file
def readFastq(filename):
    names=[]
    sequences=[]
    qualities=[]
    with open(filename) as fh:
        while True: #loop for ever
            name = fh.readline().rstrip() # read first line (name) and not save it
            seq = fh.readline().rstrip() # read second line and save to seq
            fh.readline() # read third line (+) and not save it
            qual = fh.readline().rstrip() # read fourth line (quality) and save to qual
            if len(seq) == 0: # reached end of file
                break
            names.append(name)
            sequences.append(seq)
            qualities.append(qual)
    return names, sequences, qualities

### func to create histogram for sequence qualities
def createHist(qualities):
    hist = [0] * 50 #number of bins /  max qual
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] +=1 # increment the value for that qual by one
    return hist

### func to get GC content of seqs
def findGCByPos(reads):
    gc=[0] * 100 #list for gc at each position. * 100 cause of lenght of the seqs
    totals = [0] *100 #list of totals for each pos
    
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C'or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
        for i in range(len(gc)):
            if totals[i] > 0:
                gc[i] /= float(totals[i])
    return gc 


### func to generate a naive algorithm for alignment
def naive(p, t):
    occurences=[]
    for i in range(len(t)-len(p) + 1): #loop over alignments. len(t) - len(p) + 1 ensures that it will not try to match after t ends
        match=True
        for j in range(len(p)): # loop over characters
            if t[i+j] != p[j]: # compare characters
                match=False # mismatch, reject alignment
                break
        if match:
            occurences.append(i) # record all chars matches
    return occurences

### func to generate random reads            
def generateReads(genome, numReads, readLen):
    reads= []
    for _ in range(numReads):
        start = random.randint(0, len(genome)-readLen) - 1
        reads.append(genome[start: start+readLen])
    return reads

### func to find the number of reads that match the genome
def Nmatches(reads, genome):
    numMatched=0
    n=0 # count of total number of reads processed
    for r in reads:
        r = r[:30]
        matches = naive(r, genome)
        matches.extend(naive(reverseComplement(r), genome)) # .extend adds the results of another list
        n += 1
        if len(matches) > 0: #read aligned at least in one place
            numMatched +=1
    print(' %d / %d reads matched the genome!' % (numMatched, n))
    return numMatched

### naive alignment with reverse complement
__author__ = "Ricardo Belmonte-Lopes"
def naive_with_rc(p, t):
    occurrences=[]
    f = naive(p, t) # run naive 
    r = naive(reverseComplement(p), t) #run naive with rc of p
    occurrences.extend(f) #add foward occurrences from naive
    occurrences.extend(r) #add reverse occurrences from naive with rc of p
    return sorted(occurrences)

### naive algorithm allowing two mismatches
__author__ = "Ricardo Belmonte-Lopes"
def naive_2mm(p, t):
    occurences=[]
    for i in range(len(t)-len(p) + 1): #loop over alignments. len(t) - len(p) + 1 ensures that it will not try to match after t ends
        match=True
        matches=0
        for j in range(len(p)): # loop over characters
            if t[i+j] == p[j]: # compare characters
                matches += 1 #add number of matches to matches
        if matches >= (len(p)-2): # if matches is equal or larger to the len of  p - 2 mismatches
            occurences.append(i) # record all chars matches
    return occurences

#get quals for all sequences
__author__ = "Ricardo Belmonte-Lopes"
def qualBySeq(qual):
    seq_qual_list=[]
    for qu in qual:
        seq_q=[]
        for phred in qu:
            q = phred33ToQ(phred)
            seq_q.append(q)
        seq_qual_list.append(seq_q)
    #tranform qualities in a np array and transpose to get vals of all seqs for each position   
    qual_arr= np.array(seq_qual_list)
    qual_by_pos=qual_arr.T


    #generates empty list to store mode of quals by position to make a histogram
    hist_qual=[]
    # loop over vals of all seqs for each position to get the mode for each one
    for pos in qual_by_pos:
        val, cont = stats.mode(pos)
        val_int = int(val)
        hist_qual.append(val_int)


        
    return seq_qual_list, hist_qual, qual_by_pos 

#plt.bar(range(0, len(hist_qual)), hist_qual)
#plt.show()


#Boyer-Moore learn from character comparisons to skip pointless alignments. Try alignments from left-to-right and character comparisons from right-to-left
#Bad-character rule: upon mismatch, skip alignments until mismatch becomes macht or P moves past mismatched character
#Good suffix rule: if t = substring matched by inner loop, skip until there are no mismatches between P and t or P moves past t
# Between Bad-char rule and Good suffix, uses the one that skips more

#implementing Boyer-Moore: lots of code for classes and functions in the notebook

#!/usr/bin/env python

"""bm_preproc.py: Boyer-Moore preprocessing."""

__author__ = "Ben Langmead"

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


class TestBoyerMoorePreproc(unittest.TestCase):

    def test_z_1(self):
        s = 'abb'
        #    -00
        z = z_array(s)
        self.assertEqual([3, 0, 0], z)

    def test_z_2(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_z_3(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_n_1(self):
        s = 'abb'
        #    01-
        n = n_array(s)
        self.assertEqual([0, 1, 3], n)

    def test_n_2(self):
        s = 'abracadabra'
        #    1004010100-
        n = n_array(s)
        self.assertEqual([1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11], n)

    def test_n_3(self):
        s = 'abababab'
        #    0204060-
        n = n_array(s)
        self.assertEqual([0, 2, 0, 4, 0, 6, 0, 8], n)

    def test_big_l_prime_1(self):
        s = 'abb'
        #    001
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 2], big_l_prime)

    def test_big_l_prime_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)

    def test_small_l_prime_1(self):
        s = 'abracadabra'
        # N  1004010100-
        # l'           1
        # l'        4
        # l' 44444444111
        small_l_prime = small_l_prime_array(n_array(s))
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

    def test_good_suffix_match_mismatch_1(self):
        p = 'GGTAGGT'
        big_l_prime, big_l, small_l_prime = good_suffix_table(p)
        self.assertEqual([0, 0, 0, 0, 3, 0, 0], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 3, 3, 3], big_l)
        self.assertEqual([7, 3, 3, 3, 3, 0, 0], small_l_prime)
        self.assertEqual(0, good_suffix_mismatch(6, big_l_prime, small_l_prime))
        self.assertEqual(0, good_suffix_mismatch(6, big_l, small_l_prime))
        #  t:      xT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(5, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(5, big_l, small_l_prime))
        #  t:     xGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(4, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(4, big_l, small_l_prime))
        #  t:    xGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(3, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(3, big_l, small_l_prime))
        #  t:   xAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(2, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(2, big_l, small_l_prime))
        #  t:  xTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(1, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(1, big_l, small_l_prime))
        #  t: xGTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(0, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(0, big_l, small_l_prime))

    def test_good_suffix_table_1(self):
        s = 'abb'
        #    001
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 2], big_l_prime)
        self.assertEqual([0, 0, 2], big_l)
        self.assertEqual([3, 0, 0], small_l_prime)

    def test_good_suffix_table_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        # l' -4444444111
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8], big_l)
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

if __name__ == '__main__':
    unittest.main()
""" end öf preprocessing"""


def boyer_moore(p, p_bm, t):
    i = 0 # index
    occurrences = [] # list for matches
    while i < len(t) - len(p) +1: # loop through all positions in t were p could start
        shift = 1 # indicates the amount to move along t
        mismatched = False
        for j in range(len(p) -1, -1, -1): # last -1 is to loop backwards
            if p[j] != t[i+j]: # mismatch
                skip_bc = p_bm.bad_character_rule(j, t[i+j]) # number of alignments skiped by bad character rule
                skip_gs = p_bm.good_suffix_rule(j) # number of alignments skiped by good suffix rule
                shift = max(shift, skip_bc, skip_gs) # find the maximum skip from the three variables
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip() # number of alignments to skip if match
            shift = max(shift, skip_gs)

        i += shift #update position by shift

    return occurrences
#Boyer-Moore preprocess p
# algorithm that preprocess T is called offline, otherwise algorithm is online. Multiple sequence alignment to a genome is offline

# offline preprocessessing = ordering by index or grouping by some characteristc
# generate a index for each substring from a given size (k-mer index). a substring of lenght 5 would be 5-mer index. Its a data structure called multimap
#Query index of T with substring of same lenght from p. if is in the index its a hit. Then match T and p and verificates if rest of p is match (full match) or mismatch
# the index is searched by binary search. First compares to the middle of the index. The query can be alphabetically greater or smaller (ie TGG > GTG). If larger ignore first half, if smaller ignore second. Repeat on the selected bisection. Repeat until  match. ~log2(n) bisections per query

# binary search in python
import bisect
bisect.bisect_left(a, x) # leftmost offset were x can be inserted into a to maintain order

#Hash tables: hash function h maps 3-mer to buckets

#implementing a k-mer index

class Index(object):
    def __init__(self, t, k): # t = text, k = k-mer lenght
        '''Create index from all substrings of size "lenght" '''
        self.k = k # set the k to the one passed
        self.index=[] #creates an empty list for the index
        for i in range(len(t) - k +1): # set a loop to get all k-mer and its leftmost offset
            self.index.append((t[i:i+k], i)) # append the k-mer and the index as a tuple
        self.index.sort()
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1)) # as all indexes are greater than -1, adding it ensures that it return always the first one
        hits=[] #hits of the kmer
        while i < len(self.index):
            if self.index[i][0] != kmer: # if index != kmer break cause is not in the list. the [0] indicates the kmer
                break
            hits.append(self.index[i][1]) # the one indicates the index position
            i += 1

        return hits
    


#to query the index

def queryIndex(p, t, index):
    k = index.k #lenght of k
    offsets = [] #list for the matches
    for i in index.query(p): # returns a list where the first k bases of p match the first k bases of t
        if p[k:] == t[i+k:i+len(p)]: # verification
            offsets.append(i)
    return offsets


#cut odd or even entries from index speeds the binary search. Or every third (if k =3).
#subsequence of S: string of chars also occurring in S in the same order
#seq = 'AACCGGTT'
#seq[0]+seq[1]+seq[5]+seq[7] > AAGT # subsequence
#seq.find('AAGT') > -1 #not a substring
#using subsequences increases the specificity of the filter provided by the index
  
def approximate_match(p, t, n): #n= mismatches
    seqment_lenght =  int(round(len(p) / (n+1)))
    all_matches = set()
    for i in range(n+1):
        start = i * segment_lenght
        stop = min((i+1)* segment_lenght, len(p)) # min to not allow to be larger than p
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t) # returns a list where that substring of p matched t

        #validation with no more than n mismatches
        for m in matches:
            if m < start or m-start+len(p) > len(t): # in either case p runs past the ends of t
                continue # skip the rest of this loop
            mismatches = 0
            for j in range(0, start): # test parts of p before the matched segment
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)): # test parts of p after the matched segment
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start) # m - start to get the position of the begining of p
    return list(all_matches)


### if I remember right, all of the funcs from __author__='Ben Langmead' to here were from the course lectures
                
############ modified naive function to include counts of compared alignments and characters
__author__ = "Ricardo Belmonte-Lopes"        
def naive_with_counts(p, t):
    occurences=[]
    alignments=0
    chars=0
    for i in range(len(t)-len(p) + 1): #loop over alignments. len(t) - len(p) + 1 ensures that it will not try to match after t ends
        match=True
        alignments += 1 # each try is a new alignment
        chars += 1 # for each alignment there is at least one char
        for j in range(len(p)): # loop over characters
            if t[i+j] != p[j]: # compare characters
                
                match=False
                #
                # mismatch, reject alignment
                break
        chars += j  # add j number of tries when break    
        
        # chars += 1 #gives char =1
        if match:
            occurences.append(i) # record all chars matches
    return occurences, alignments, chars

#######################################################


############ modified boyer moore function to include counts of compared alignments and characters
__author__ = "Ricardo Belmonte-Lopes"
def boyer_moore_with_counts(p, p_bm, t):
    i = 0 # index
    occurrences = [] # list for matches
    alignments=0
    chars= 0

    while i < len(t) - len(p) +1: # loop through all positions in t were p could start
        alignments += 1 # add alignment 
        shift = 1 # indicates the amount to move along t
        mismatched = False
        for j in range(len(p) -1, -1, -1): # last -1 is to loop backwards
            chars +=1 # add char comparison
            if p[j] != t[i+j]: # mismatch
                skip_bc = p_bm.bad_character_rule(j, t[i+j]) # number of alignments skiped by bad character rule
                skip_gs = p_bm.good_suffix_rule(j) # number of alignments skiped by good suffix rule
                shift = max(shift, skip_bc, skip_gs) # find the maximum skip from the three variables
                mismatched = True
                break
        
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip() # number of alignments to skip if match
            shift = max(shift, skip_gs)

        i += shift #update position by shift

    return occurrences, alignments, chars

############################


############## get all kmers of k lenght for p ################
__author__ = "Ricardo Belmonte-Lopes"
def getKmer(p, k):
    p_kmer_list=[]
    for x in range(len(p)-k+1):
        kmer = p[x:x+k]
        s = kmer, x
        p_kmer_list.append(s)
    return p_kmer_list


####### get all indexes for P up to n mismatches from a dictIndex object (definition below)
__author__ = "Ricardo Belmonte-Lopes"
def query_dictIndex(slice_d, p_kmer_list, n):
    new_dict={}
    for key, pos in slice_d.items():
        for j in range(len(p_kmer_list)):
            mismatches=0
            for i in range(len(key)):
                if p_kmer_list[j][0][i] != key[i]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                new_dict[key]={}
                
                new_dict[key]['start']=0
                new_dict[key]['start']= new_dict[key]['start']+p_kmer_list[j][1]
                new_dict[key]['pos']=[]
                if type(pos) == list:
                    new_dict[key]['pos'].extend(pos)
                else:
                    new_dict[key]['pos'].append(pos)
    return new_dict

######### modified Index class using dicts and allowing mismatches to generate query base
__author__ = "Ricardo Belmonte-Lopes"        
class dictIndex(object):
    def __init__(self, t, k): # t = text, k = k-mer lenght, n = mismatches allowed
        '''Create index from all substrings of size "lenght" '''
        self.k = k # set the k to the one passed
        self.dict={} #generates an empty dict to store the main index
        for i in range(len(t) - k +1): # set a loop to get all k-mer and its leftmost offset
            if t[i:i+k] not in self.dict.keys(): # if a kmer is not in dict, add it as key, and add i to that key
                self.dict[t[i:i+k]]=[]
                self.dict[t[i:i+k]].append(i)
            else: 
                self.dict[t[i:i+k]].append(i) # if key already in dict, just add i to key
    def getKmer(p, n):
        p_kmer_list=[]
        for x in range(len(p)-n+1):
            kmer = p[x:x+n]
            s = kmer, x
            p_kmer_list.append(s)
        return p_kmer_list

    def query(self, p, n): #method to query the index. need a kmer list generated by getKmer(). Allow n mismatches
        k=self.k
        p_kmer_list = getKmer(p, k)
        new_dict={} # generates a new dict to store the results
        for key, pos in self.dict.items(): # loop over the keys and items in the main index dict
            for j in range(len(p_kmer_list)): #for each of the kmers in the list
                mismatches=0 #start a mismatches object
                for i in range(len(key)): # create a range of the length of the key to slice the key and the kmer
                    if p_kmer_list[j][0][i] != key[i]: # if i letter of key != of i letter of p[j][0] (0 indicates the string), add 1 to mismatches
                        mismatches += 1
                        if mismatches > n: # break
                            break #changed from pass to break
                if mismatches <= n: 
                    new_dict[key]={} # create a new key = key in the new dict
                    new_dict[key]['start']=0 # creates a subkey start (to store the start position in relation to p) in the previous key, and initiate it to 0.
                    new_dict[key]['start']= new_dict[key]['start']+p_kmer_list[j][1] # add the start position related to p in the start subkey
                    new_dict[key]['pos']=[] #creates a subkey pos in the previous created key, it will store all pos where that key appers in the object used to create the index
                    if type(pos) == list: #if the pos are of type list, use extend
                        new_dict[key]['pos'].extend(pos)
                    else: 
                        new_dict[key]['pos'].append(pos) # if pos are of type int, use append
        
            
        return new_dict
          
                        
    # Final query for the dict index with mismatches                
    def queryIndex(self, p, t, n):
        st= d.datetime.now()
        k = self.k #lenght of k
        offsets = [] #list for the matches
        q_dict = self.query(p, n) #generates a query dict using the query method from the index, using a p kmers list of p generated by getKmer(), and  allowing n mismatches
        total_index_hits = 0 # to store number of total index hits
        for key, items in q_dict.items(): # iterates over query dict
            
            for i in range(len(q_dict[key]['pos'])): # for each position in pos for a given kmer key
                if q_dict[key]['pos'][i] < len(p): # if the position is smaller than p, break
                    break
                else:
                    x = q_dict[key]['pos'][i] - q_dict[key]['start'] # get the correct start position accounting for the positions of the different kmers of p
                    y = q_dict[key]['pos'][i] + (len(p) - q_dict[key]['start']) ## get the correct stop position accounting for the positions of the different kmers of p
                    sub = t[x:y] # slice t with the calculated start and stop, generating a subtring
                    #mismatches = 0
                    if sub == p: # if the substring == to p then add its start position to offsets and total index hits
                        offsets.append(x)
    
                    else:
                        mismatches=0 # initiates a mismatches object, which will reset at each new position in a kmer key
                        for j in range(len(sub)): # generate a range of len(sub) to slice p and sub
                            if p[j] != sub[j]: # if p != sub, add 1 to mismatches
                                mismatches += 1
                                if mismatches > n: # if mismatches greater than n break
                                    break
                if mismatches <= n : #if after all j possibilities mismatches still <= n, add sub start to offsets
                    offsets.append(x)

            
            offsets = sorted(list(set(offsets))) #remove duplicates of offsets using set() and return it to list format
        # to register the total number of index hits
        for key in q_dict.keys():
            total_index_hits += len(q_dict[key]['pos'])
            
        return offsets, total_index_hits, print(d.datetime.now()-st).total_seconds()
                   





#############################modified to method query return index hits up to n mismatches
__author__ = "Ricardo Belmonte-Lopes"
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    

    def genSubseq(self, p):
        k = self.k
        ival=self.ival
        span = 1 + ival * (k - 1)
        subseq=[]
        for x in range(len(p)-span+1):
            pos_p = x
            slice_p = p[x:x+span:ival]
            comb = slice_p, pos_p
            subseq.append(comb)
        return subseq

    def query(self, p, t, n):
        """ Return index hits for first subseq of p """
        subseq = self.genSubseq(p)
        hits = []

        for s in range(len(subseq)):

            for i in range(len(self.index)):
                mismatches=0
                for j in range(len(self.index[i][0])):
                    
                    if self.index[i][0][j] != subseq[s][0][j]:
                        mismatches+=1
                        if mismatches > n:
                            break
                if mismatches <= n:
                    
                    pos = self.index[i][1] -s
                    sub_t = t[pos:pos+(len(p)-s)]
                    for o in range(len(sub_t)):
                        if sub_t[o] != p[o]:
                            mismatches += 1
                            if mismatches > n:
                                break
                    if mismatches <= n:
                        hits.append(pos)
        
        hits_final=list(set(hits))
            

        return hits


########################################


#generate subseq index for p
__author__ = "Ricardo Belmonte-Lopes"
def genSubseq_Index(p, k,ival):
    k = k # number of chars
    ival = ival #spaceing between chars
    span = 1 + ival * (k - 1)
    subseq=[]
    for x in range(len(p)-span+1):
        pos_p = x
        slice_p = p[x:x+span:ival]
        comb = slice_p, pos_p
        subseq.append(comb)
    return subseq


  

        
##### function to answer HW2Q5 ########
__author__ = "Ricardo Belmonte-Lopes"
def index_approximate_match(p, t, index, n): #n= mismatches
    segment_lenght =  int(round(len(p) / (n+1)))
    all_matches = set()
    matches_count=0 # object to count total index hits
    for i in range(n+1):
        start = i * segment_lenght
        stop = min((i+1)* segment_lenght, len(p)) # min to not allow to be larger than    
        matches = index.query(p[start:stop]) # returns a list where that substring of p matched the index
        
        #validation with no more than n mismatches
        for m in matches:
            matches_count += 1 # increase number of total index hits
            if m < start or m-start+len(p) > len(t): # in either case p runs past the ends of t
                continue # skip the rest of this loop
            mismatches = 0
            for j in range(0, start): # test parts of p before the matched segment
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(stop, len(p)): # test parts of p after the matched segment
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start) # m - start to get the position of the begining of p
                #matches_count += 1
    return list(all_matches), matches_count



###########original######

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

##### func to answer HW2Q6######
__author__ = "Ricardo Belmonte-Lopes"    
def Subindex_approximate_match(p, t, subindex, n): #n= mismatches
    segment_lenght =  int(round(len(p) / (n+1)))
    all_matches = set()
    matches_count=0
    ival = subindex.ival
    span = subindex.span
    
    for i in range(n+1):

        start = i * segment_lenght
        stop = min((i+1)* segment_lenght, len(p)) # min to not allow to be larger than p
        sub_p_ = p[i:span+i:ival]
        matches =[]
        # just to count full hits on the index
        for w in range(len(subindex.index)):
            if subindex.index[w][0] == sub_p_:
                matches_count +=1
            
        # to acctually get all occurrences up to 2 mismatches
        sub_p = (p[start:stop:ival])
        for j in range(len(subindex.index)):
            if sub_p in subindex.index[j][0]:
                matches.append(subindex.index[j][1])
       
        
        #validation with no more than n mismatches
        for m in matches: 
            mismatches = 0
            for l in range(len(p[start:stop])):  # test if index hit matched segment has mismatches
                if p[start:stop][l] != t[m:(m+len(p[start:stop]))][l]:
                    mismatches += 1
                    if mismatches > n:
                        break

            for j in range(0, start): # test parts of p before the matched segment
                if start == 0: # if start in first pos of p break
                    break
                else:
                    if p[j] != t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break
            for j in range(stop, len(p)): # test parts of p after the matched segment
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m-start) # m - start to get the position of the begining of p
                #matches_count += 1 #counts only matches that passed the mismatches validation
    return list(all_matches), matches_count

#Hamming distance: minimum number of substitutions need to turn one seq in the other (if same lenght)
#Levenshtein distance: minimum number of edits (substitutions, indels)


### func to find the max distance between seqs using a Naive Hamming algorithm
def naiveHamming(p, t, maxDistance):
    occurences=[]
    for i in range(len(t)-len(p) + 1): #loop over alignments. len(t) - len(p) + 1 ensures that it will not try to match after t ends
        nmm=0
        match=True
        for j in range(len(p)): # loop over characters
            if t[i+j] != p[j]: # compare characters
                nmm +=1 # mismatch
                if nmm > maxDistance:
                    break # exceeded max hamming dist
        if nmm <= maxDistance:
            occurences.append(i) # approximate match
    return occurences

### A recursive func to find edit distances between two seqs
def edDistRecursive(a,b):
    st=d.datetime.now()
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    delt = 1 if a[-1] != b[-1] else 0
    return min(edDistRecursive(a[:-1], b[:-1]) + delt, edDistRecursive(a, b[:-1]) +1, edDistRecursive(a[:-1], b) +1), print(d.datetime.now()-st).total_seconds()

### func to find edit distance
def edDistance(x, y):
    D = [] # empty 'matrix'
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1)) #creates the matrix and fill with zeroes
    for i in range(len(x)+1):
        D[i][0] = i #fills first line
    for i in range(len(y)+1):
        D[0][i] = i # fills first row
    #fills the rest of the matrix row  by row, column by column
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] +1 # edit dist from cell to left plus 1
            distVer = D[i-1][j] + 1 #edit dist from cell above plus 1
            if x[i-1] == y[j-1]: #if char of above and left cells are the same
                distDiag = D[i-1][j-1] #value of diagonal is the same from left cell of the diagonal
            else:
                distDiag = D[i-1][j-1] +1 # if the chars dont match add one

            D[i][j] = min(distHor, distVer, distDiag) #minimize the distance for the cell

    return D[-1][-1] , D # returns bottom right value, which is the minimum edit distance
        
#%%time #clock a funct in Jupyter

##edit distance and approximate matching matrix: T is the row, P is the column
# first row is initialized with all zeroes (represents offsets, all equally probable) and first column with 0 to len(P)
## min edit distance will be in final row. Look for origin of the value using the quadrants(hor, ver, diag) [traceback]. Vertical move indicates insertion between P and T. 
# using a penalty matrix for transitions (2), transversion (4) and gaps or indels (8) can be accomplished by changing the added term of each D (hor, ver, diag)
# global alignment : for the whole seq / local alignmnt: finds the most similar pair of substrings from X and Y
# for local alignment we search for max distance, and instead of a penalty matrix is used an scoring matrix (matches are positive, differences are negative). Traceback using the local alignment matrix stops if reaches a zero


#global alignment

def globalAlignment(x, y):
    alphabet = ['A', 'T', 'C', 'G']
    score = [[0, 4, 2, 4, 8], # five by five, one for each char and one for skip
             [4, 0, 4, 2, 8], 
             [2, 4, 0, 4, 8], 
             [4, 2, 4, 0, 8], 
             [8, 8, 8, 8, 8]]
    D = [] # empty 'matrix'
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1)) #creates the matrix and fill with zeroes
        
    for i in range(1, len(x)+1): #range to left the left/top zero
        D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1] #find first char in the row and the value for skiping
        
    for i in range(len(y)+1):
        D[0][i] = D[0][i-1] + score[-1][alphabet.index(y[i-1])] # fills first row
    #fills the rest of the matrix row  by row, column by column
        
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]#  # pound if skip char in y
            distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1] ##pound if skip char in x
            if x[i-1] == y[j-1]: #if char of above and left cells are the same
                distDiag = D[i-1][j-1] #value of diagonal is the same from left cell of the diagonal
            else:
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])] # score for mismatch

            D[i][j] = min(distHor, distVer, distDiag) #minimize the distance for the cell

    return D[-1][-1] # returns bottom right value, which is the minimum edit distance

#assembly: compares prefix of one seq with suffix of another. 1st law: if a suffix of read A is similar to a prefix of read B then A and B might overlap in the genome. Second law: more coverge lead to more and longer overlaps


#find overlap between two strings

def overlap(a, b, min_overlap=3):
    start= 0 #index

    while True: # while there is overlap of mininum size
        start = a.find(b[:min_overlap], start) # find first overlap of minimum size and set the new start
        if start == -1: # no overlap
            return 0
        if b.startswith(a[start:]): # verification of the overlap
            return len(a) - start # lenght of the overlap
        start += 1 # if none of the conditions meet increment start and reenter loop

# extended overlap func to create an overlap map

def overlap(a, b, min_overlap):
    """"Return lenght of longest suffix of a matching a prefix of b that is at least min overlap lenght long. If no overlap exists return 0"""
    start= 0 #index, start at leftmost pos

    while True: # while there is overlap of mininum size
        start = a.find(b[:min_overlap], start) # find first overlap of minimum size and set the new start
        if start == -1: # no overlap, no more occurences to right
            return 0
        if b.startswith(a[start:]): # verification if found an overlap
            return len(a) - start # lenght of the overlap
        start += 1 # if none of the conditions meet increment start and reenter loop

def naive_overlap_map(reads, k):
    olaps={}
    for a, b in permutations(reads, 2): # for each permutations of two seqs from reads
        olen= overlap(a, b, min_overlap=k) #run overlap func
        if olen > 0: #if there is overlap
            olaps[(a,b)] = olen # add tuple a,b as key to dict and overlap lenght as value
    return olaps
            
######for HW3

### func to find edit distance between seqs
__author__ = "Ricardo Belmonte-Lopes"
def edDistance_apr(x, y):
    D = [] # empty 'matrix'
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1)) #creates the matrix and fill with zeroes
    for i in range(len(x)+1):
        D[i][0] = i #fills first line
    for i in range(len(y)+1):
        D[0][i] = 0 # fills first column
    #fills the rest of the matrix row  by row, column by column
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] +1 # edit dist from cell to left plus 1
            distVer = D[i-1][j] + 1 #edit dist from cell above plus 1
            if x[i-1] == y[j-1]: #if char of above and left cells are the same
                distDiag = D[i-1][j-1] #value of diagonal is the same from left cell of the diagonal
            else:
                distDiag = D[i-1][j-1] +1 # if the chars dont match add one

            D[i][j] = min(distHor, distVer, distDiag) #minimize the distance for the cell

    return min(D[-1]) # returns bottom right value, which is the minimum edit distance
    

# get all kmers from a list of reads
__author__ = "Ricardo Belmonte-Lopes"
def kmerDict(p_kmer_set):
    kmerdict=dict.fromkeys(p_kmer_set)
    for key in kmerdict.keys():
        kmerdict[key]=set()
    return kmerdict

### get overlap between seqs

def overlap(a, b, min_overlap):
    """"Return lenght of longest suffix of a matching a prefix of b that is at least min overlap lenght long. If no overlap exists return 0"""
    start= 0 #index, start at leftmost pos

    while True: # while there is overlap of mininum size
        start = a.find(b[:min_overlap], start) # find first overlap of minimum size and set the new start
        if start == -1: # no overlap, no more occurences to right
            return 0
        if b.startswith(a[start:]): # verification if found an overlap
            return len(a) - start # lenght of the overlap
        start += 1 # if none of the conditions meet increment start and reenter loop

### map overlaps with naive algorithm
__author__ = "Ricardo Belmonte-Lopes"        
def naive_overlap_map(reads, k):
    olaps={}
    for a, b in permutations(reads, 2): # for each permutations of two seqs from reads
        olen= overlap(a, b, min_overlap=k) #run overlap func
        if olen > 0: #if there is overlap
            olaps[(a,b)] = olen # add tuple a,b as key to dict and overlap lenght as value
    return olaps

### read FastQ sequences modified
__author__ = "Ricardo Belmonte-Lopes"
def readFastq(filename):
    names=[]
    sequences=[]
    qualities=[]
    with open(filename) as fh:
        while True: #loop for ever
            name = fh.readline().rstrip() # read first line (name) and not save it
            seq = fh.readline().rstrip() # read second line and save to seq
            fh.readline() # read third line (+) and not save it
            qual = fh.readline().rstrip() # read fourth line (quality) and save to qual
            if len(seq) == 0: # reached end of file
                break
            names.append(name)
            sequences.append(seq)
            qualities.append(qual)
    return names, sequences, qualities

### find all overlaps between reads
__author__ = "Ricardo Belmonte-Lopes"
def ovlp_all(reads, k):
    st= d.datetime.now() # start clocking
    kmerdict={} #generate an empty dict for the kmers
    oldict={} # generate an empty dict for the overlaps
    suffix_count=set() # generate an empty set to add reads with suffix overlap
    #to get all kmers from all seqs
    for i in range(len(reads)): # for a given sequence in the file
        for j in range(len(reads[i])-k+1): # for the lenght of that sequence - lenght of kmer+1
            kmer = reads[i][j:j+k] # slice the kmers
            if kmer not in kmerdict.keys(): # if the kmer not in the kmerdict add it
                kmerdict[kmer]=set()
                kmerdict[kmer].add(reads[i]) #add reads that match that kmer
            else:
                kmerdict[kmer].add(reads[i]) #add reads that match that kmer if already in the dict
    #reduce the dict removing all entries with only one match (seq match to itself only)
    new_dict = {key:val for key, val in kmerdict.items() if len(val) > 1} # filter kmers with only one match

    for read in reads: #loop the seqs in the file
        sf=read[-k:] #slice the last k positions, the suffix
        
        if sf not in new_dict.keys(): # if the suffix is not in the reduced dict
            pass # do nothing
        else:
            l = new_dict[sf] #slice the new dict with the suffix and put result in a list
            for item in l: # for each item of that list
                if item != read: # if the item is different from the read (no need to compare the read with itself)
                    ovlen=overlap(read, item, k) # calculates the overlap using the overlap func
                    if ovlen >= k: # if overlap equal or larger than minimal
                        oldict[(read, item)] = ovlen #add entry to overlap dict
                        suffix_count.add(read) # add the read to the suffix set
    return oldict, suffix_count, print((d.datetime.now()-st).total_seconds())          

            
#import site
#site.addsitedir('/the/path')
#import sys    



####4th week


## generate a list of all shortest superstrings. Modified from scs function
__author__ = "Ricardo Belmonte-Lopes"
def scs_list(ss): ##find shortest common superstring
    shortest_sup = None #initalize shortest superstring object
    s_list=[]
    shortest_sub=set()
    s = ss.copy()
    for s_perm in itertools.permutations(s): #create all permutations for strings in ss
        sup = s_perm[0]
        for i in range(len(s)-1):
            ov_len = overlap(s_perm[i], s_perm[i+1], 1) #find overlap between two strings
            sup += s_perm[i+1][ov_len:] #append suffix
            
        if shortest_sup is None or len(sup) <= len(shortest_sup): # if shortest sup was not actualized (still none), actualize, other compare if the length of the sup is smaller than the already recorded
            shortest_sup = sup
            s_list.append(shortest_sup)
    for i in range(len(s_list)):
        if len(s_list[i]) == len(s_list[-1]):
            shortest_sub.add(s_list[i])
    return  shortest_sub #, shortest_sup #uncomment to return just one shortest superstring

################################################
#find maximum overlap

def pic_max_ov(reads,k):
    reada, readb = None, None
    best_ov_len = 0
    for a, b in itertools.permutations(reads, 2):
        ov_len = overlap(a, b, k)
        if ov_len > best_ov_len:
            reada, readb = a, b
            best_ov_len = ov_len
    return reada, readb, best_ov_len

####################################################################

#de novo assembly for HW4 Q3-4. Started from greedy shortest common superstring algorithm but ended all diferent, not used pic_max_overlap func
__author__ = "Ricardo Belmonte-Lopes"
def greedy_scs(seqs, k):
    st= d.datetime.now() # start clocking
    reads=seqs.copy()

    while len(reads) > 1:
        kmerdict={} #generate an empty dict for the kmers
        oldict={} # generate an empty dict for the overlaps
        for i in range(len(reads)): # for a given sequence in the file
            for j in range(len(reads[i])-k+1): # for the lenght of that sequence - lenght of kmer+1
                kmer = reads[i][j:j+k] # slice the kmers
                if kmer not in kmerdict.keys(): # if the kmer not in the kmerdict add it
                    kmerdict[kmer]=set()
                    kmerdict[kmer].add(reads[i]) #add reads that match that kmer
                else:
                    kmerdict[kmer].add(reads[i]) #add reads that match that kmer if already in the dict
        #reduce the dict removing all entries with only one match (seq match to itself only)
        new_dict = {key:val for key, val in kmerdict.items() if len(val) > 1} # filter kmers with only one match

        for read in reads: #loop the seqs in the file
            max_ov_len=[]
            sf=read[-k:] #slice the last k positions, the suffix
            
            if sf not in new_dict.keys(): # if the suffix is not in the reduced dict
                pass # do nothing
            else:
                l = new_dict[sf] #slice the new dict with the suffix and put result in a list
                for item in l: # for each item of that list
                    if item != read: # if the item is different from the read (no need to compare the read with itself)
                        ovlen=overlap(read, item, k) # calculates the overlap using the overlap func
                        if ovlen >= k: # if overlap equal or larger than minimal
                            if read not in oldict.keys():
                                oldict[read]= {}
                                oldict[read]['overlap']=[]
                                oldict[read]['overlap'].append((item, ovlen)) #add entry to overlap dict
                            else:
                                oldict[read]['overlap'].append((item, ovlen))


        for key in oldict.keys():
            max_ov=[] # create an empty list to add the maximum overlap values for comparison
            oldict[key]['max overlap']=[] #creates an entry in the dict to add the maximum overlap seq and value
            for i in range(len(oldict[key]['overlap'])): #for all overlaps of a given key
                max_ov.append(oldict[key]['overlap'][i][1]) #append the overlap value to the list
                m = max(max_ov) #find the maximum overlap value
                idx_m = max_ov.index(m) #get the index of the maximum overlap value
                read_x = oldict[key]['overlap'][idx_m][0] #slice the key for the seq with that index
            oldict[key]['max overlap'].append((read_x, m)) #append the seq and the overlap value to the dict
            
        klist=[] #creates an empty list to store less than ideal overlaps
        for key in oldict.keys():
            comp = oldict[key]['max overlap'][0][0] # get the seq of maximum overlap for a given key
            if comp in oldict.keys(): # if the seq above is in the dict keys, compare to see if it have a larger overlap with other sequence. If not add the seq to the klist
                if oldict[key]['max overlap'][0][1] >= oldict[comp]['max overlap'][0][1]:
                    klist.append(comp)



        for i in range(len(klist)):
            ld = list(oldict.keys())
            if klist[i] in ld: # if a sequence in the klist is in the dict keys, delete from dict
                del oldict[klist[i]]

                    
        for key in oldict.keys():
            cut = oldict[key]['max overlap'][0][1] #get value of overlap for a key
            reads.append(key+oldict[key]['max overlap'][0][0][cut:]) #append to the reads the seq combined with the one of maximum overlap 
            if key in reads:
                reads.remove(key) #remove original seq from the reads 
            if oldict[key]['max overlap'][0][0] in reads:
                reads.remove(oldict[key]['max overlap'][0][0]) # remove the overlaped seq from the reads
            
                
    return reads[0], print((d.datetime.now()-st).total_seconds())
                

                      
###################
__author__ = "Ricardo Belmonte-Lopes"
def reduce(seqs, k):
    reads=seqs.copy()
    kmerdict={} #generate an empty dict for the kmers
    oldict={} # generate an empty dict for the overlaps
    for i in range(len(reads)): # for a given sequence in the file
        for j in range(len(reads[i])-k+1): # for the lenght of that sequence - lenght of kmer+1
            kmer = reads[i][j:j+k] # slice the kmers
            if kmer not in kmerdict.keys(): # if the kmer not in the kmerdict add it
                kmerdict[kmer]=set()
                kmerdict[kmer].add(reads[i]) #add reads that match that kmer
            else:
                kmerdict[kmer].add(reads[i]) #add reads that match that kmer if already in the dict
    #reduce the dict removing all entries with only one match (seq match to itself only)
    new_dict = {key:val for key, val in kmerdict.items() if len(val) > 1} # filter kmers with only one match

    for read in reads: #loop the seqs in the file
        max_ov_len=[]
        sf=read[-k:] #slice the last k positions, the suffix
        
        if sf not in new_dict.keys(): # if the suffix is not in the reduced dict
            pass # do nothing
        else:
            l = new_dict[sf] #slice the new dict with the suffix and put result in a list
            for item in l: # for each item of that list
                if item != read: # if the item is different from the read (no need to compare the read with itself)
                    ovlen=overlap(read, item, k) # calculates the overlap using the overlap func
                    if ovlen >= k: # if overlap equal or larger than minimal
                        if read not in oldict.keys():
                            oldict[read]= {}
                            oldict[read]['overlap']=[]
                            oldict[read]['overlap'].append((item, ovlen)) #add entry to overlap dict
                        else:
                            oldict[read]['overlap'].append((item, ovlen))


    for key in oldict.keys():
        max_ov=[] # create an empty list to add the maximum overlap values for comparison
        oldict[key]['max overlap']=[] #creates an entry in the dict to add the maximum overlap seq and value
        for i in range(len(oldict[key]['overlap'])): #for all overlaps of a given key
            max_ov.append(oldict[key]['overlap'][i][1]) #append the overlap value to the list
            m = max(max_ov) #find the maximum overlap value
            idx_m = max_ov.index(m) #get the index of the maximum overlap value
            read_x = oldict[key]['overlap'][idx_m][0] #slice the key for the seq with that index
        oldict[key]['max overlap'].append((read_x, m)) #append the seq and the overlap value to the dict
        
    klist=[] #creates an empty list to store less than ideal overlaps
    for key in oldict.keys():
        comp = oldict[key]['max overlap'][0][0] # get the seq of maximum overlap for a given key
        if comp in oldict.keys(): # if the seq above is in the dict keys, compare to see if it have a larger overlap with other sequence. If not add the seq to the klist
            if oldict[key]['max overlap'][0][1] >= oldict[comp]['max overlap'][0][1]:
                klist.append(comp)



    for i in range(len(klist)):
        ld = list(oldict.keys())
        if klist[i] in ld: # if a sequence in the klist is in the dict keys, delete from dict
            del oldict[klist[i]]

                
    for key in oldict.keys():
        cut = oldict[key]['max overlap'][0][1] #get value of overlap for a key
        reads.append(key+oldict[key]['max overlap'][0][0][cut:]) #append to the reads the seq combined with the one of maximum overlap 
        if key in reads:
            reads.remove(key) #remove original seq from the reads 
        if oldict[key]['max overlap'][0][0] in reads:
            reads.remove(oldict[key]['max overlap'][0][0]) # remove the overlaped seq from the reads
        
                
    return reads

                        
                        
  
########################################


    

#de bruijn graph method

def de_bruijn(st, k):
    edges=[]
    nodes=set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1], st[i+1:i+k])) #get right and left kmers
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

import graphviz

def plot_de_bruijn(st, k):
    nodes, edges = de_bruijn(st, k)
    dot_str = 'digraph "De Bruijn: graph" {\n'
    for node in nodes:
        dot_str += ' %s (label="%s"); \n'
    for src, dst in edges:
        dot_str += '%s -> %s; \n'
    return dot_str + '}\n'


    
