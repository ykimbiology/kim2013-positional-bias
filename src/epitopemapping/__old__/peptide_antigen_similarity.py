#! /usr/bin/python

import re

'''
    [20111011] This script appears redundant and confusing. Remove it if not being used.
    Appears to be bruteforce approach to peptide:protein mapping.
    
 [] Remove redundant code:
   When is a peptide present in an antigen?
   When mapping epitopes of OrganismSource to antigens of OrganismTarget, we have to define when we can map these peptides.

   + Name sould be changed to peptide_antigen_mapping.py
   + May need to have its own module for future expansion.
   '''

class PeptideMatch(object):
    def __init__(self):
        self.is_match = False
        self.match_result_list = [] # [fraction_sim, seq_peptide, index_a, index_b]

    def get_match_state(self):
        if len(self.match_result_list) > 0: self.is_match=True
        return self.is_match

def pepmatch(seq_peptide, seq_antigen,fraction_sim_cutoff=0.80):
    '''Searches for a similar peptide in a given antigen. Returns true if there is a match.'''
    pm = PeptideMatch()
    pm.match_result_list = map_epitope_to_antigen_fraction_sim(seq_peptide, seq_antigen, fraction_sim_cutoff = fraction_sim_cutoff) # [fraction_sim, seq_peptide, index_a, index_b]
    return pm

def map_peptides(dic_peptides, antigen_list, dic_homology, fraction_sim_cutoff=0.80, debug=False):
    '''Simplify the process of mapping. Separate assigning weights. Make it work with other programs.'''
    #fraction_sim_cutoff=0.80
    dic_content_stats = {} # Also include information for unmapped_orfs

    peptide_id_list = dic_peptides.keys()

    if debug==True: print 'len(dic_peptides)', len(dic_peptides)
    if debug==True: print dic_peptides

    for (i, record) in enumerate(antigen_list): # For each antigen, return a list of peptides that are found.
        id_antigen  = record.id.split('|')[1].strip()  # Assumes the format: gi|###|XXXXXX
        #description = record.description.split('|')[4]
        description = '|'.join(record.description.split('|')[2:])
        #description = record.description
        seq_antigen = str(record.seq)
        match_list_temp = None

        if dic_homology != None:
            match_list_temp = [pepmatch(dic_peptides[epitope_id]['linear_peptide_seq'],seq_antigen,fraction_sim_cutoff=fraction_sim_cutoff) for epitope_id in peptide_id_list if dic_homology[epitope_id].__contains__(id_antigen)==True] # Loops over all epitopes; attemsp matching only those whose source_antigen are homologous.
        else:
            match_list_temp = [pepmatch(dic_peptides[epitope_id]['linear_peptide_seq'],seq_antigen,fraction_sim_cutoff=fraction_sim_cutoff) for epitope_id in peptide_id_list] # Loops over all epitopes; attemsp matching only those whose source_antigen are homologous.
        match_list = [m for m in match_list_temp if m.get_match_state()==True]
        if debug==True: print 'map_peptides, len(match_list)', i, len(match_list)
        if len(match_list) > 0:
            if debug==True: print '\tmap_peptides', '\t', i, '\t', len(match_list), '\t', id_antigen, '\t', record.description
            dic_content_stats[id_antigen] = {'id_antigen':id_antigen, 'match_list':match_list, 'description':description, 'num_peptides':len(match_list), 'mapped':True}
        else:
            dic_content_stats[id_antigen] = {'id_antigen':id_antigen, 'match_list':None, 'description':description, 'num_peptides':0, 'mapped':False}
    return dic_content_stats


def map_peptides_b(dic_peptides, antigen_list, dic_homology, fraction_sim_cutoff=0.80, debug=False):
    '''This version bypasses object creation;Simplify the process of mapping. Separate assigning weights. Make it work with other programs.'''
    #fraction_sim_cutoff=0.80
    dic_content_stats = {} # Also include information for unmapped_orfs

    peptide_id_list = dic_peptides.keys()

    if debug==True: print 'len(dic_peptides)', len(dic_peptides)
    if debug==True: print dic_peptides

    for (i, record) in enumerate(antigen_list): # For each antigen, return a list of peptides that are found.
        id_antigen  = record.id.split('|')[1].strip()  # Assumes the format: gi|###|XXXXXX
        #description = record.description.split('|')[4]
        description = '|'.join(record.description.split('|')[2:])
        #description = record.description
        seq_antigen = str(record.seq)
        match_list_temp = None

        if dic_homology != None:
            match_list_temp = [map_peptide(dic_peptides[epitope_id]['linear_peptide_seq'],seq_antigen,fraction_sim_cutoff=fraction_sim_cutoff) for epitope_id in peptide_id_list if dic_homology[epitope_id].__contains__(id_antigen)==True] # Loops over all epitopes; attemsp matching only those whose source_antigen are homologous.
        else:
            match_list_temp = [map_peptide(dic_peptides[epitope_id]['linear_peptide_seq'],seq_antigen,fraction_sim_cutoff=fraction_sim_cutoff) for epitope_id in peptide_id_list] # Loops over all epitopes; attemsp matching only those whose source_antigen are homologous.
        match_list = [m for m in match_list_temp if len(m)>0]
        if debug==True: print 'map_peptides, len(match_list)', i, len(match_list)
        if len(match_list) > 0:
            if debug==True: print '\tmap_peptides', '\t', i, '\t', len(match_list), '\t', id_antigen, '\t', record.description
            dic_content_stats[id_antigen] = {'id_antigen':id_antigen, 'match_list':match_list, 'description':description, 'num_peptides':len(match_list), 'mapped':True}
        else:
            dic_content_stats[id_antigen] = {'id_antigen':id_antigen, 'match_list':None, 'description':description, 'num_peptides':0, 'mapped':False}
    return dic_content_stats






def get_peptides(seq_antigen, peptide_len):
    peptide_list = []
    num_peptides = len(seq_antigen) - peptide_len + 1
    for i in range(num_peptides):
        peptide = seq_antigen[i:i + peptide_len]
        peptide_list.append(peptide)
    return peptide_list

def get_hamming_fraction(peptide_a, peptide_b):
    '''Fraction of residues shared between the two peptides.
       One way to speed up is to use dot product of vectors? and sum 1s? Not sure how much faster?'''
    count_shared = 0
    for i in range(len(peptide_a)):
        res_a = peptide_a[i]
        res_b = peptide_b[i]
        if (res_a == res_b):
            count_shared = count_shared + 1
    fraction_shared = float(count_shared)/len(peptide_a)
    return fraction_shared

def map_epitope_to_antigen(seq_peptide, seq_antigen):
    '''For now, find an exact match between the peptide and the antigen.
       '''
    is_match = False
    index_a = -1
    index_b = -1

    match_result_list = []
    seq_peptide = seq_peptide.upper()
    seq_antigen = seq_antigen.upper()

    p = re.compile(seq_peptide)
    mobj = re.search(p, seq_antigen)
    if (mobj != None):
        is_match = True
        index_a = mobj.start()
        index_b = mobj.end()
    #result = (index_a, index_b)
    result = (index_a, index_b, seq_peptide)
    match_result_list.append(result)
    return (is_match, match_result_list)


def map_epitope_to_antigen_fraction_sim(seq_peptide, seq_antigen, fraction_sim_cutoff = 0.80):
    '''Given a peptide and an antigen, returns a list of matched peptides and their corresponding positions.'''
    #print 'debug map_epitope_to_antigen_fraction_sim', seq_peptide, seq_antigen
    match_result_list = []
    if ((seq_antigen != None) and (seq_peptide != None)):
        peptide_list = get_peptides(seq_antigen, len(seq_peptide))
        fraction_sim_list = [get_hamming_fraction(peptide, seq_peptide) for peptide in peptide_list]
        for pos in range(len(peptide_list)):
            peptide = peptide_list[pos]  # This peptide is from seq_anatigen, NOT the query peptide.
            fraction_sim = fraction_sim_list[pos]
            if (fraction_sim >= fraction_sim_cutoff):
                index_a = pos
                index_b = pos + len(peptide)
                #result = [fraction_sim, peptide, index_a, index_b]
                result = [fraction_sim, seq_peptide, index_a, index_b]  # Q: Should query_peptide be used here? or matched_peptide?
                match_result_list.append(result)
        if len(match_result_list) == 0: match_result_list = []
    return match_result_list



