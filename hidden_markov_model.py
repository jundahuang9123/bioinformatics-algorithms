#!/usr/bin/env python3

"""
Author: Junda Huang 910203370050

Description: this is a script to implement Hidden Markov Model
"""
# Import statements
from sys import argv
import random

# Function definitions

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 } 
    
def parse_fasta(filename):
    """
    Parse sequences from fasta file into dictionary

    input:
        filename: (string) filename of the fast file
    output:
        seqdict: (dictionary) dictionary with sequences as value
            title as key. e.g. seqdict = {>HBA_HUMAN: VFA--HAGEY}
    """
    seqdict = {}
    with open(filename, 'r') as fasta:
        for line in fasta:
            a = str(line).strip()
            if a[0] == '>':
                seq = ''
                key = a
                seqdict[key] = []
            else:
                seq += a
                seqdict[key] = seq
    return seqdict      

def match_state_info(seqdict, threshold):
    """
    find out the number of match states and their positions 
        with give set of sequences

    input:
        seqdict: (dictionary) dictionary with sequences as value
            title as key. e.g. seqdict = {>HBA_HUMAN: VFA--HAGEY}
        threshould: (float) a number where how many percentage of amindacid 
            present count as a match state
    output:
        nmatcches: (int) match states number
        match_pos: (list)  of int of match states positions in the sequneces
    """
    nmatches = 0
    match_pos = []
    seqlist = [seq for key, seq in seqdict.items()]
    for i in range(len(seqlist[0])):
        count = 0
        for j in range(len(seqlist)):
            if seqlist[j][i] != '-':
                count +=1
        if count/len(seqlist) >= threshold:
            nmatches += 1
            match_pos.append(i)
    return nmatches, match_pos

def hmm_state_check(seqdict, nmatches, match_pos, pa):
    """
    check each position of sequences for their state

    input:
        seqdict: (dictionary) dictionary with sequences as value
            title as key. e.g. seqdict = {>HBA_HUMAN: VFA--HAGEY}
        nmatches: (int) match states number
        match_pos: (list)  of int of match states positions in the sequneces
        pa: (dictionary) Background amino acid probabilities
    output:
        states_dict: (dictionary) of states and pos in dictionary as value,
            sequences title as key, position as subkey 
            e.g. states_dict = {>HBA_HUMAN: {'begin': M, 1: I, 2: D...}
    """
    states_dict = {}
    for key, seq in seqdict.items():
        states_dict[key] = {}
        if seq[0] in pa.keys() and 0 not in match_pos:
            states_dict[key]['begin'] = 'I'
        else:
            states_dict[key]['begin'] = 'M'
        for i in range(1, len(seq)+1):
            if i < len(seq):
                if seq[i] in pa.keys() and i in match_pos:
                    states_dict[key][i-1] = 'M'
                elif seq[i] in pa.keys() and i not in match_pos:
                    states_dict[key][i-1] = 'I'
                elif seq[i] not in pa.keys() and i in match_pos:
                    states_dict[key][i-1] = 'D'
            elif i == len(seq):
                states_dict[key]['end'] = states_dict[key][i-2]
    return states_dict

def transitions(states_dict):
    """
    count the number of trainsitions between states for all sequences
        and for each sequences 

    input:
        states_dict: (dictionary) of states and pos in dictionary as value,
            sequences title as key, position as subkey 
            e.g. states_dict = {>HBA_HUMAN: {'begin': M, 1: I, 2: D...}
    output:
        n_transitions: (int) total number of transitions
        transition_list: (list) of int of transitions in sequences
    """
    n_transitions = 0
    transition_list = []
    for states in states_dict:
        tr = 0
        value = list(states_dict[states].values())
        for i, pos in enumerate(value):
            if i == 0 :
                continue
            elif pos != value[i-1]:
                n_transitions += 1
                tr += 1
        transition_list.append(tr)    
    return n_transitions, transition_list

def transition_probabilities(match_pos, states_dict):
    """
    calculate transitions probabilities at each position

    input: 
        match_pos: (list)  of int of match states positions in the sequneces  
        states_dict: (dictionary) of states and pos in dictionary as value,
            sequences title as key, position as subkey 
            e.g. states_dict = {>HBA_HUMAN: {'begin': M, 1: I, 2: D...}
    output:
        tran_p_dict: (dictionary) state as key of dictionaries
            position as subkey probabilities and posistions in tuple as value, 
            e.g. tran_p_dict = {I: {0: [(I, 0.5), (M, 0.5), (D, 0)]},
                                   {1: [(M, 0), (D, 0), (M, 0)]}...}
    """
    tran_p_dict = {} # construting HMM dictionary of dictionaries
    tran_p_dict['M'] = {}; tran_p_dict['I'] = {}; tran_p_dict['D'] = {}
    for key in tran_p_dict:
        for i in range(len(match_pos)+1):
            if i == 0:
                tran_p_dict[key]['begin'] = []
            else:
                tran_p_dict[key][i-1] = []
    values = list(states_dict.values())
    keyset = []
    for value in values: 
        keyset += [list(value.keys())]
    match_pos_new = [x-1 for x in match_pos]
    for j in range(len(match_pos)+1):
        dict_trans = {'MM': 0, 'MI': 0, 'MD': 0, 'IM': 0, \
            'II': 0, 'ID': 0, 'DM': 0, 'DI': 0, 'DD': 0}   
    # initiallise probabilities dataset in another dictionary for calculation
        for k, l in enumerate(keyset):
            if j == 0: # check begin-0 states change
                key = 'M' + values[k][l[j]]
                if values[k][l[j]] == 'I':
                    dict_trans['I'+values[k][l[j+1]]] += 1
            else:
                if j+1 < len(l):
                    key = values[k][l[j]] + values[k][l[j+1]]
                    # convert states into dataset key
            dict_trans[key] += 1
        sum_m = dict_trans['MM']+dict_trans['MI']+dict_trans['MD']
        sum_i = dict_trans['IM']+dict_trans['II']+dict_trans['ID']
        sum_d = dict_trans['DM']+dict_trans['DI']+dict_trans['DD']
        # sum up for M, I and D
        sumlist = [sum_m, sum_i, sum_d]
        for m, state in enumerate(tran_p_dict):
        # calculation and store data accordingly
            if j == 0:
                subkey = 'begin'
            else:
                subkey = j - 1
            if sumlist[m] == 0:
                tran_p_dict[state][subkey] += [('M', 0), ('I', 0), ('D', 0)]
            else:
                keyone = ''.join(map(str, state))
                tran_p_dict[state][subkey] += \
                    [('M', dict_trans[keyone+'M']/sumlist[m]), \
                        ('I', dict_trans[keyone+'I']/sumlist[m]), \
                            ('D', dict_trans[keyone+'D']/sumlist[m])] 
    return tran_p_dict

def match_emissions(seqdict, match_pos):
    """
    find all residues that are emitted in match states

    input:
        seqdict: (dictionary) dictionary with sequences as value
            title as key. e.g. seqdict = {>HBA_HUMAN: VFA--HAGEY}
        match_pos: (list)  of int of match states positions in the sequneces
    outpus:
        emissions_m: (list) of list of emitted residues in position order
    """
    emissions_m = []
    for i, seq in enumerate(seqdict.values()):
        for j, pos in enumerate(match_pos):
            if i == 0 :
                emissions_m.append([])
            if seq[pos] in pa.keys():
                emissions_m[j].append(seq[pos])
    return emissions_m

def emission_probabilities_m(emissions_m, pa, pseudocounts = False):
    """
    calculate all emission probabilities in match states

    input:
        emissions_m: (list) of list of emitted residues in position order
        pseudocounts: (boolean) set critereia if pseudocounts required = True
        pa: (dictionary) Background amino acid probabilities 
    output:
        emission_p_m: (dictionary) with match state positions as key, 
            AAs and probabilities in tuple as value; e.g.
            emission_p_m = {0: (A, 0.75), (N, 0.125), (M, 0.125)}
    """
    emission_p_m = {}
    for i, aalist in enumerate(emissions_m):
        l = len(aalist)
        aas = set(aalist) # find unique amino acids present in values
        emission_p_m[i] = []
        if pseudocounts == True:
            for aa in pa.keys:
                emission_p_m[i] += [(aa, (aalist.count(aa)+1)/(l+20))]
        elif pseudocounts == False:
            for aa in aas:
                emission_p_m[i] += [(aa, aalist.count(aa)/l)]
    return emission_p_m

def pssm_format(emission_p_m):
    """
    print a matrix of emission probabilities in PSSM format

    input:
        emission_p_m: (dictionary) with match state positions as key, 
            AAs and probabilities in tuple as value; e.g.
            emission_p_m = {0: (A, 0.75), (N, 0.125), (M, 0.125)}
    output:
        print out e.g. 
            	    A       R	    N       D	    C	  
                1 A	2.507	-1.007	-1.001	-1.526	-0.387	
                2 A	-1.379	-2.418	-3.295	-3.531	-0.261	
                3 A	1.286	-1.389	3.736	-0.498	-0.573
    """
    header = ' \t \t' + '\t'.join(aa for aa in pa.keys())
    matrix_main = [header]
    for key, value in emission_p_m.items():
        linestart = str(key) + ' A\t'
        line = linestart
        for aa in pa.keys():
            for i in range(len(value)):
                prob = 0
                if value[i][0] == aa:
                    prob = value[i][1]
                    break
            line = line + '\t' + '{:.3f}'.format(prob)    
        matrix_main.append(line)
    for matrix_line in matrix_main:
        print(matrix_line)
    return 

def sequence_generator(emission_p_m, tran_p_dict, pa, rep):
    """
    generating possible sequences based on hmm profile

    input:
        emission_p_m: (dictionary) with match state positions as key, 
            AAs and probabilities in tuple as value; e.g.
            emission_p_m = {0: (A, 0.75), (N, 0.125), (M, 0.125)}
        tran_p_dict: (dictionary) state as key of dictionaries
            position as subkey probabilities and posistions in tuple as value, 
            e.g. tran_p_dict = {I: {0: [(I, 0.5), (M, 0.5), (D, 0)]},
                                   {1: [(M, 0), (D, 0), (M, 0)]}...}
        pa: (dictionary) Background amino acid probabilities
        rep: (int) amount of sequences required
    output:
        sequences_gen: (list) of generated sequences
        sequences_states: (list) of states of sequences of their postions
    """
    sequences_gen = []
    sequences_states = []
    for h in range(rep):
        seq_len_count = 0
        sequence = ''
        sequence_state = ''
        m = 0
        state_next = ''.join(map(str, random.choices(['I', 'M'], \
            [tran_p_dict['M']['begin'][1][1], \
                tran_p_dict['M']['begin'][0][1]], cum_weights = None, k = 1)))
        # find first residue of the sequences
        if state_next == 'I':
            aa = random.choices(list(pa.keys()), weights = \
                    [k for k in pa.values()], cum_weights = None, \
                        k = 1)
            sequence = sequence+ ' '.join(map(str, aa))
        elif state_next == 'M' and tran_p_dict['M']['begin'][0][1] > 0:
            sequence = sequence + '-' 
        else:
            aa = random.choices(emission_p_m[m], \
                weights = [l[1] for l in emission_p_m[m]], \
                    cum_weights = None, k = 1)[0][0]
            sequence = sequence + (' '.join(map(str, aa)))
            m += 1
            sequence_state = sequence_state + state_next
        while state_next != 'end':
            weights = []
            for n in range(len(tran_p_dict)):
                weights += [tran_p_dict[state_next][m][n][1]]
                # assign probabilities from dictionary for choices 
            state_next = ''.join(map(str, random.choices(list(\
                tran_p_dict.keys()), weights, cum_weights = None, k = 1)))
            # choose next state 
            if state_next == 'D':
                sequence = sequence + '-'
                m += 1
            elif state_next == 'I':
                aa = random.choices(list(pa.keys()), weights = \
                    [k for k in pa.values()], cum_weights = None, k = 1)
                sequence = sequence + ' '.join(map(str, aa))
            elif state_next == 'M':
                aa = random.choices(emission_p_m[m], \
                    weights = [l[1] for l in emission_p_m[m]], \
                        cum_weights = None, k = 1)[0][0]
                sequence = sequence + (' '.join(map(str, aa)))
                m += 1
            sequence_state = sequence_state + state_next
            if m == len(tran_p_dict['M'])-1:
                sequences_gen.append(sequence)
                sequences_states.append(sequence_state)
                state_next = 'end'
    return sequences_gen, sequences_states

if __name__ == "__main__":

    # implement main code here
    #infile = 'test.fasta'
    infile = 'test.fasta'
    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
    
    # Question 1:
    seqdict = parse_fasta(infile) 
    # Parse sequences data 
    nmatches, match_pos = match_state_info(seqdict, threshold = 0.5)
    # Get match states information
    print('Answer to question 1:\n{} match states needed.'.format(nmatches))
    
    # Question 2:
    states_dict = hmm_state_check(seqdict, nmatches, match_pos, pa)
    # check the states of each residues in each sequences 
    tran_p_dict = transition_probabilities(match_pos, states_dict)
    # calculate each state transition probabilities in each position
    print('\nAnswer to question 2:\nTransition probability are as follow:')
    for i, v in tran_p_dict.items():
        for j in v.items():
            print(i, j)
    emissions_m = match_emissions(seqdict, match_pos)
    # find out all match state emissions
    emission_p_m = \
        emission_probabilities_m(emissions_m, pa, pseudocounts = False)
    # calculate match state emission probabilities per position
    print('\nThe emission probabilities are as follow:\n') 
    for k, v in emission_p_m.items():
        print(k, v)
    
    # Question 3:
    print('\nAsnwer to question 3:')
    pssm_format(emission_p_m)
    # convert emission probabilities into PSSM format and print
    
    # Question 4:
    sequences_gen, sequences_states = \
        sequence_generator(emission_p_m, tran_p_dict, pa, rep = 10)
    print('\nAnswer to question 4:\nThe generated sequences are as follow:')
    for i in sequences_gen:
        print(i)
    
    # Question 5:
    infile_large = 'test_large.fasta'
    # 1:
    seqdict_large = parse_fasta(infile_large) 
    nmatches_large, match_pos_large = \
        match_state_info(seqdict_large, threshold = 0.5)
    print('Answer to question 1:\n{} match states needed.'\
        .format(nmatches_large))
    # 2:
    states_dict_large = hmm_state_check\
        (seqdict_large, nmatches_large, match_pos_large, pa)
    tran_p_dict_large = transition_probabilities\
        (match_pos_large, states_dict_large)
    print('\nAnswer to question 2:\nTransition probability are as follow:')
    for i_large, v_large in tran_p_dict_large.items():
        for j_large in v_large.items():
            print(i_large, j_large)
    emissions_m_large = match_emissions(seqdict_large, match_pos_large)
    emission_p_m_large = \
        emission_probabilities_m(emissions_m_large, pa, pseudocounts = False)
    print('\nThe emission probabilities are as follow:\n', emission_p_m_large)
    # 3:
    print('\nAsnwer to question 3:')
    pssm_format(emission_p_m_large)
    # 4:
    sequences_gen_large, sequences_states_large = \
        sequence_generator(emission_p_m_large,\
             tran_p_dict_large, pa, rep = 10)
    print('\nAnswer to question 4:\nThe generated sequences are as follow:')
    for i_large in sequences_gen_large:
        print(i_large)