#!/usr/bin/env python3
"""
Author: Junda Huang
Student number: 910203370050
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table 
for a query sequence
- write a function to find the best hit
- write a fasta parser to read in the Arabidopsis data

"""
# import statements
from sys import argv
from urllib.request import urlopen
import re

# implement your functions here
def parse_fasta(filename):
    """
    Parse DNA sequences from fasta file into a list

    input:
        filename: (string) utf coded fasta file
    return:
        seqlist: (list) a list of DNA sequences
        seqscount: (int) amount of sequences in the fasta
        totallen: (int) total length of the sequence(s)
    """
    seqlist = []
    filelist = []
    seqstartlist = []
    totallen = 0
    for line in filename: # decode fasta from URL
        decoded_line = line.decode("utf-8")
        a = str(decoded_line)
        filelist.append(a) # make into a list of strings per line
    # make the whole list into a long string, take out "'", ",", and " "
    seqliststr = str(filelist).replace("'", "")\
        .replace(",", "").replace(" ", "")
    for b in range(len(seqliststr)):
        if seqliststr[b] == ">": # find sequences start positions
            seqstartlist += [b]
    seqscount = len(seqstartlist) # count amount of sequences
    for c in range(seqscount):
        if c+1 != seqscount: # except last sequences
            for d in range(seqstartlist[c], seqstartlist[c+1]):
                if seqliststr[d:d+2] == "\\n": 
                    # find the first linebreaker before the next sequence start 
                    # then start record sequence
                    seqlist += [seqliststr[d+2:seqstartlist[c+1]]\
                        .replace("\\n", "")]
                    break
        else:
            for e in range(seqstartlist[c], len(seqliststr)):
                if seqliststr[e:e+2] == "\\n":
                    # find the first linebreaker then start record sequence
                    seqlist += [seqliststr[e+2:]\
                        .replace("\\n", "")]
                    break
    for f in seqlist: # calculate total sequence length
        totallen += len(f)
    return seqlist, seqscount, totallen

def hash(seqd, k):
    """
    Creating hash tabels with given k-tuple and given set of sequences D

    input:
        seqd: (list) a list of DNA sequences
        k: (int) the length of the k-tuple
    return:
        hashdict: (dictionary) dictionary of hash table, the k-tuples as keys
                    postions as tuple as (seq1, pos1)
    """
    hashdict = {}
    for h, seq in enumerate(seqd): # find k-tuples that existed in all sequences
        for i in range(0, len(seq), k):
            if seq[i:i+k] not in hashdict.keys():
                hashdict[seq[i:i+k]] = []
                hashdict[seq[i:i+k]] += [(h+1, i+1)]
            elif 'N' in seq[i:i+k]: # ignore the 'N'
                continue
            else:
                hashdict[seq[i:i+k]].append((h+1, i+1))
    return hashdict

def seq_search(seqq, hashdict, k, ends = False):
    """
    Compare given query sequence to the hash table to get amount of hits

    input: 
        seqq: (string) a sequence to be queried
        hashdict: (dictionary) the hash table to be queried on
        k: (int) the length of the k-tuple
        ends: (boolean) true or false, if want to extend the beginning and end 
    return:
        hits: (list) a list of tuples with the hits information
            [(index1, shift1, offset1), (index2, shift2, offset2), ...]
    """
    subseqq = [] 
    # create a list of seqq subsequences in length of k-tuple
    for i in range(len(seqq)-k+1):
        subseqq += [seqq[i:i+k]]
    subseqq += [seqq[len(seqq)-k+1:]]
    hits = [] 
    # create a list of hits
    if ends == False:
        for j, t in enumerate(subseqq):
            if t in hashdict.keys():
                for pos in hashdict[t]:
                    hits += [(pos[0], pos[1]-j, pos[1])]
    elif ends == True: 
    # for extension of begining and end in case of less than k match
        for j, t in enumerate(subseqq):
            if j == 0 or j == len(subseqq)-1:
                if '\t'.join(t) in hashdict.keys():
                # begin and end with substring matches
                    for key in hashdict:
                        if t in key:
                            for pos in hashdict[key]:
                                hits += [(pos[0], pos[1]-j, pos[1])]     
            else:
                if t in hashdict.keys():
                    for pos in hashdict[t]:
                        hits += [(pos[0], pos[1]-j, pos[1])]
    hits.sort()
    return hits

def run(hits):
    """
    Find out the longest "run" in the list of hits

    input: 
        hits: (list) the hits list of tuples
    return:
        runmaxlist: (list) the longest run(s) in a list 
        runmax: (int) longest run length
    """
    hitdict = {} 
    # Covert hit list of tuples into a dictionary of dictionary
    runmax = 0
    runmaxlist = []
    for hit in hits:
        if hit[0] not in hitdict.keys(): 
            hitdict[hit[0]] = {} # the primary key is the index
            hitdict[hit[0]][hit[1]] = 1 
            # sub key is the shift; value is the repetition of the shift
        elif hit[1] not in hitdict[hit[0]].keys():
            hitdict[hit[0]][hit[1]] = 1
        else:
            hitdict[hit[0]][hit[1]] += 1
    for key in hitdict: # find longest run number
        for subkey in hitdict[key]:
            if hitdict[key][subkey] >= runmax:
                runmax = hitdict[key][subkey]
    for key in hitdict: # find longest run position
        for subkey in hitdict[key]:
            if hitdict[key][subkey] == runmax:
                for hit in hits:
                    if key == hit[0] and subkey == hit[1]:
                        runmaxlist += [(key, subkey, hit[2])]
    return runmaxlist, runmax

def align_graph(runmaxlist, runmax, seqd, seqq, k):
    """
    Create sequences alignment graphs

    input:
        runmaxlist: (list) of hits which have longest "run"
        runmax: (int) of the length of the longest "run"
        seqd: (list) sequences were queried on
        seqq: (string)query sequence 
        k: (int) k-tuple length
    return:
        alignmentlist: (list) list of (lists of) listsstrings 
                        with aligned sequences and alignment symbols
    """
    seqquery = ''
    seqalign = ''
    alignmentlist = []
    for i in range(0, len(runmaxlist), runmax):
        alignmentlist += [[]]
        seqnum = runmaxlist[i][0]-1
        # find the target sequences no. in seqd
        posstart = runmaxlist[i][2]
        # find the starting point of alignment
        posend = runmaxlist[i+runmax-1][2]+k-1
        # find the ending point of alignment
        seqquery = ' '*(posstart-1) + seqq + \
            ' '*(len(seqd[seqnum])-posend) 
        # contruct query seq with spaces
        if seqq == seqd[seqnum]: # contruct alignment with spaces, | and .
            seqalign = ' '*(posstart-1) + '|'*len(seqq) + \
                ' '*(len(seqd[seqnum])-posend)
        else:
            for j in range(len(seqq)):
                if seqq[j] != seqd[seqnum][j+posstart-1]:
                    seqalign += '.'
                else:
                    seqalign += '|'
            seqalign = ' '*(posstart-1) + seqalign + \
                ' '*(len(seqd[seqnum])-posend)
        alignmentlist += [[seqquery, seqalign, seqd[seqnum]]] 
        # add 3 alignment line into a list of a list
    return alignmentlist

def rev_com(seqs):
    """
    Create a list of reverse complementary sequences

    input:
        seqs: (list) DNA seqeunces
    output:
        revseqs: (list) reverse complementary DNA sequences
    """
    nucs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    # create dictionary and a pattern for subsitution
    nucs = dict((re.escape(ori), rep) for ori, rep in nucs.items())
    pattern = re.compile('|'.join(nucs.keys())) 
    revseqs = []
    for seq in seqs:
        # subsitutions
        newseq = pattern.sub(lambda m: nucs[re.escape(m.group(0))], seq)
        reseq = ''.join(reversed(newseq)) #reverse
        revseqs += [reseq]
    return revseqs

if __name__ == "__main__":

    # the code below should produce the results necessary to answer 
    # the questions. In other words, if we run your code, we should 
    # see the data that you used to answer the questions
    
    query = 'TGCAACAT'
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1,s2,s3]

    url1 = 'http://www.bioinformatics.nl/courses/BIF-31306/TAIR10.fasta'
    url2 = 'http://www.bioinformatics.nl/courses/BIF-31306/query.fasta'
    file_d = urlopen(url1) # obtain TAIR10.fasta file
    file_query = urlopen(url2) # obtain query.fasta file

    # Question 1: print hash table 
    print("Answer to Q1:")
    hashdict = hash(seqs, k = 2)
    hashdictlist =[]
    for i in hashdict:
        hashdictlist = hashdict[i]
        print (i, ' '.join(map(str, hashdictlist)))

    # Question 2: print hit result
    hits = seq_search(query, hashdict, k = 2)
    print("\nAnswer to Q2:\n", len(hits), hits[0], hits[-1])

    # Question 3: print alignment graph
    runmaxlist, runmax = run(hits)
    print("\nAnswer to Q3:\n", runmax, ' '.join(map(str, runmaxlist)))
    alignmentlist = align_graph(runmaxlist, runmax, seqs, query, k = 2)
    for i in alignmentlist:
        for j in i:
            print(j)

    # Question 4: Parsing Arabidopsis thaliana genome sequences 
    seq_arabi, seqs_count_arabi, len_arabi = parse_fasta(file_d)
    print("\nAnswer to Q4:\nThere are {1} sequences in the genome file.\n\
        The total length of the genome sequence is {0} bp."\
            .format(len_arabi, seqs_count_arabi))

    # Question 5: build hash tabel for Arabidopis chromosomes
    hash_arabi = hash(seq_arabi, k = 13)
    print("\nAnswer to Q5:\nThere are {} keys in my hash table/dictionary."\
        .format(len(hash_arabi)))

    # Question 6: parse query sequences and find hits in Arabidopis hash table
    queries, seqs_count_query, len_query = parse_fasta(file_query)
    runresult = {}
    for queryarabi in queries: 
    # This should have been done earlier. Realised later that 
    # the seq_search function only takes a single query sequence at once
        hitsarabi = seq_search(queryarabi, hash_arabi, k = 13)
        runmaxlistara, runmaxara = run(hitsarabi)
        runresult[runmaxara] = runmaxlistara
    result = max(runresult)
    print("\nAnswer to Q6:\n Maximum hits is {}.\nThe hits are {}"\
        .format(result, ' '.join(map(str, runresult[result]))))

    # Optional Question: finding hits of reverse complement genome
    rearabi = rev_com(queries) # create reverse complement query sequences
    runresultre = {}
    for queryarab in rearabi:
    # same issue as in Question 6
        hitsarabire = seq_search(queryarab, hash_arabi, k = 13, ends = True)
        runmaxlistarare, runmaxarare = run(hitsarabire)
        runresultre[runmaxarare] = runmaxlistarare
    resultre = max(runresultre)
    print("\nAnswer to Optional:\n Maximum hits is {}.\nThe hits are {}"\
        .format(resultre, ' '.join(map(str, runresultre[resultre]))))