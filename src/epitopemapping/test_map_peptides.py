#! /usr/bin/python

import time
import os
import cPickle
import sys

#TODO: Need to fix these hardcoded paths.
sys.path.append('/home/yohan/workspace/immunomebrowser/core/src/immunomebrowser.main/immunomebrowser')
sys.path.append('/home/yohan/workspace/immunomebrowser/core/src/immunomebrowser.epitopemapping/immunomebrowser')

"""
Make sure peptide:protein mapping has been done correctly.
Compare fast methods against bruteforce approach, which is a lot slower.
"""
from Bio import SeqIO

from mappeptides import map_peptides_selector

DATASET_NAME   = 'vacv'
PATH_MAPPING   = '/home/yohan/data_large/immunomebrowser/experiments/compare_alignment_algorithm/20120730_test_usearch/mapping/'+DATASET_NAME

#== VACV:
PATH_DATA_TEST = '/home/yohan/data_large/immunomebrowser/experiments/test_data/sequences/vacv'
FNAME_PEPTIDE_LIST = os.path.join(PATH_DATA_TEST, 'peptide_list.2k.fasta')
FNAME_PROTEIN_LIST = os.path.join(PATH_DATA_TEST, 'sequences.10245.fasta')

#== HCV:
#PATH_DATA_TEST = '/home/yohan/data_large/immunomebrowser/experiments/test_data/sequences/hcv'
#FNAME_PEPTIDE_LIST = os.path.join(PATH_DATA_TEST, 'peptide_list.2k.fasta')
#FNAME_PROTEIN_LIST = os.path.join(PATH_DATA_TEST, 'sequences.11103.fasta')

def get_d_peptide_sequence():
    d={}
    for seq_record in SeqIO.parse(FNAME_PEPTIDE_LIST, "fasta"):
        id_peptide = seq_record.id.replace('id|','')
        peptide = str(seq_record.seq)
        d[id_peptide] = peptide
        #print seq_record.id
        #print str(seq_record.seq)
        #print len(seq_record)
    return d

def generate_d_mapping(method = 'bruteforce', fname_peptide_list='peptide_list.fasta', sim_cutoff_fraction=1.0):
    '''test_map_peptides: Map peptides using a given algorithm.'''
    fname_peptide_list = FNAME_PEPTIDE_LIST
    fname_antigens_target = FNAME_PROTEIN_LIST
    
    ti = time.time()
    (d_mapping, mappingAlgo) = map_peptides_selector(fname_peptide_list, fname_antigens_target, method=method, sim_cutoff_fraction=sim_cutoff_fraction)
    tj = time.time()
    td = tj-ti
    
    fname_output = '.'.join(['d_mapping',method, str(sim_cutoff_fraction), 'cpickle'])
    fname_output = os.path.join(PATH_MAPPING, fname_output)
    f=open(fname_output,'w')
    cPickle.dump(d_mapping, f)
    f.close()
    
    print method, len(d_mapping)
    print 'TIME ', method,  ' %.2f seconds' %(td,)
    #assert False

    
def get_statistics_count_shared(id_list_a, id_list_b):
    '''Return number of entries shared between the two lists.
       And its fraction out of the larger list.
    '''
    num_a = len(id_list_a)
    num_b = len(id_list_b)
    num_max = max((num_a, num_b))
    id_list_intersection = list(set(id_list_a).intersection(set(id_list_b)))
    id_list_difference   = list(set(id_list_a).difference(set(id_list_b)))   # A list of ids in 'A' that are NOT in 'B':
    num_shared = len(id_list_intersection)
    num_missed = len(id_list_difference)
    
    num_shared_fraction = num_shared/float(num_max)
    #print 'num_shared_fraction', [num_shared_fraction, num_a, num_b, num_shared, num_missed]
    return num_shared_fraction, id_list_intersection, id_list_difference

def get_statistics_count_shared_proteins(id_list_protein_a, id_list_protein_b, output=True):
    #print 'id_list_protein_a', sorted(id_list_protein_a)[0:5]
    #print 'id_list_protein_b', sorted(id_list_protein_b)[0:5]
    (num_shared_fraction, id_list_intersection, id_list_diff) = get_statistics_count_shared(id_list_protein_a, id_list_protein_b)
    num_shared_fraction_str = "%.3f" %(num_shared_fraction,)
    if output==True:
        print '== proteins with mapped peptides =='
        print '    method_a', len(id_list_protein_a)
        print '    method_b', len(id_list_protein_b)
        print '    shared  ', len(id_list_intersection)
        print '    shared_f', num_shared_fraction_str
    row = [len(id_list_protein_a), len(id_list_protein_b), len(id_list_intersection), num_shared_fraction_str]
    return row
    
    
def get_statistics_count_shared_peptides(d_mapping_a, d_mapping_b, output=True):
    ''' '''
    d_mapping_peptide_a = get_d_mapping_peptides(d_mapping_a)
    d_mapping_peptide_b = get_d_mapping_peptides(d_mapping_b)
    id_list_peptide_a = d_mapping_peptide_a.keys()
    id_list_peptide_b = d_mapping_peptide_b.keys()
    
    #print 'DEBUG ', id_list_peptide_a[0:10]
    #print 'DEBUG ', id_list_peptide_b[0:10]
    
    (num_shared_fraction, id_list_intersection_peptide, id_list_diff) = get_statistics_count_shared(id_list_peptide_a, id_list_peptide_b)
    num_shared_fraction_str = "%.3f" %(num_shared_fraction,)
    #print '\t'.join(map(str,['Fraction of peptides shared (f_shared, count_a, count_b, intersect):', num_shared_fraction_str, len(id_list_peptide_a), len(id_list_peptide_b), len(id_list_intersection_peptide)]))
    
    
    '''For the intersection of peptide_list, find out whether mapped regions are the same.'''
    count_mapping_shared = 0
    for id_peptide in id_list_intersection_peptide:
        #== GOAL: Find out for each mapping, the first one is contained in the second one, rather than whether they are equal.
        if d_mapping_peptide_a[id_peptide] == d_mapping_peptide_b[id_peptide]:
            count_mapping_shared = count_mapping_shared + 1
    mapping_peptide_shared_fraction = float(count_mapping_shared)/len(id_list_intersection_peptide)
    mapping_peptide_shared_fraction_str = "%.3f" %(mapping_peptide_shared_fraction,)
    if output==True:
        print '== peptides that were mapped =='
        print '    method_a', len(id_list_peptide_a)
        print '    method_b', len(id_list_peptide_b)
        print '    shared  ', len(id_list_intersection_peptide)
        print '    shared_f', num_shared_fraction_str, '# Fraction taken with respect to the larger set of peptides.'
        print '    shared_f_mapping_region', mapping_peptide_shared_fraction_str, '#For each id_peptide, do two sets of mapped region match?' #, len(id_list_intersection_peptide), count_mapping_shared
    
    d_peptide_seq = get_d_peptide_sequence()
    output_missed_peptides = False
    if output_missed_peptides == True:
        id_list_diff = sorted(id_list_diff, key = lambda x: int(x))
        for id_peptide in id_list_diff:
            print id_peptide, d_peptide_seq[id_peptide]
        
    row = [len(id_list_peptide_a), len(id_list_peptide_b), len(id_list_intersection_peptide), num_shared_fraction_str, mapping_peptide_shared_fraction_str] 
    return row
    
    
def get_d_mapping_peptides(d_mapping):
    '''
    d[id_peptide] = a list of (sim, id_protein, pos_start, pos_end)
    
    (sim, id_peptide, pos_start, pos_end)
    (sim, id_protein, pos_start, pos_end)
    
    '''
    d_mapping_peptide = {}
    
    id_list_peptide = []
    id_list_protein = d_mapping.keys()
    for id_protein in id_list_protein:
        mapping_list = d_mapping[id_protein] 
        for (sim, id_peptide, pos_start, pos_end) in mapping_list:
            row = (sim, id_protein, pos_start, pos_end)
            d_mapping_peptide[id_peptide] = sorted(d_mapping_peptide.get(id_peptide, []) + [row])
            
    return d_mapping_peptide


def read_d_mapping(method='blast_a', sim_cutoff_fraction=1.0):
    fname = '.'.join(['d_mapping', method, str(sim_cutoff_fraction), 'cpickle'])
    fname = os.path.join(PATH_MAPPING, fname)
    d_mapping = cPickle.load(open(fname,'r'))
    return d_mapping

def count_match_with_partial_query_peptide(d_mapping):
    '''
    Q: For a given mapping results, how many partial matches are there?
    OUTPUT: total_matches, count_matches_partial_query
    '''
    d_peptide = get_d_peptide_sequence()
    count_matches_total = 0
    count_matches_w_partial_query = 0
    
    id_antigen_list = d_mapping.keys()
    for id_antigen in id_antigen_list:
        mlist = d_mapping[id_antigen]
        for m in mlist:
            count_matches_total = count_matches_total + 1
            (sim, id_peptide, pos_i, pos_j) = m
            plen = pos_j - pos_i + 1
            peptide = d_peptide[id_peptide]
            if plen == len(peptide):
                count_matches_w_partial_query = count_matches_w_partial_query + 1
    return (count_matches_total, count_matches_w_partial_query)
            

def compare_mapping_algorithms(method_a='bruteforce', method_b='blast_a', sim_cutoff_fraction=1.0):
    ''' (1) Compare two sets of id_proteins: bruteforce and blast.  What fraction is shared?
        
       Q: What is the condition that tells us if two mappings are the same?
       A: (1) For each protein, the same set of peptides are found.
          (2) For each peptide, same number of hits are found.
          (3) For each mapping, same coordinates are found.
    RETURN: 
        (1) Percentage of shared number of proteins:  num_proteins_shared/min(num_protein_total_a, num_protein_total_b)
        (2) Percentage of shared number of peptides: num_peptides_shared/max(num_peptides_total_a, num_peptides_total_b)
        (3) Of those shared peptides, fraction of peptides with the same mapping results:
        
    RESULTS:
        == bruteforce vs. blast_a
        Fraction of proteins shared    0.927272727273    165    159    153
        Fraction of peptides shared    0.907251264755    593    540    538
        mapping_peptide_shared_fraction 0.992565055762 538 534
        
        == bruteforce vs. blast_b
        Fraction of proteins shared    0.846153846154    165    195    165
        Fraction of peptides shared    0.985049833887    593    602    593
        mapping_peptide_shared_fraction 0.954468802698 593 566
        
        == bruteforce vs. usearch_ublast
        Fraction of proteins shared (f_shared, count_a, count_b, intersect):    0.981927710843    165    166    163
        Fraction of peptides shared (f_shared, count_a, count_b, intersect):    0.954468802698    593    567    566
        mapping_peptide_shared_fraction 0.994699646643 566 563

        == bruteforce vs. smithwaterman
        Fraction of proteins shared    0.570934256055    165    289    165
        Fraction of peptides shared    0.973727422003    593    609    593
        mapping_peptide_shared_fraction 0.836424957841 593 496

       '''
    
    #== Read in previously generated mapping results using diff methods:
    d_mapping_a = read_d_mapping(method=method_a, sim_cutoff_fraction=sim_cutoff_fraction)
    d_mapping_b = read_d_mapping(method=method_b, sim_cutoff_fraction=sim_cutoff_fraction)
    
    #== statistics:
    row_protein = get_statistics_count_shared_proteins(d_mapping_a.keys(), d_mapping_b.keys(), output=False)
    row_peptide = get_statistics_count_shared_peptides(d_mapping_a, d_mapping_b, output=False)
    [count_total, count_partial] = count_match_with_partial_query_peptide(d_mapping_b)
    row = row_protein+row_peptide+[count_total, count_partial]
    #print '\t'.join(map(str, row))
    return row


def get_method_list(): 
    #return ['bruteforce', 'blast_a', 'ublast',  'search_global_fulldp', 'search_local', 'search_local_cmplt', 'usearch_local']
    return ['bruteforce', 'blast_a',  'search_global_fulldp', 'search_local', 'search_local_cmplt', 'usearch_local']
    
def get_cutoff_list(): return [1.0, 0.8, 0.6, 0.4]
    



if __name__ == '__main__':
    
    #== Generating Mapping Data:
#    method_list = get_method_list()
#    cutoff_list = get_cutoff_list()
#    for cutoff in cutoff_list:
#        for method in method_list[1:]:
#        #for method in ['bruteforce']:
#            print [cutoff, method]
#            generate_d_mapping(method = method, sim_cutoff_fraction=cutoff)
    
    
    method_ref = 'bruteforce'
    method_list = get_method_list()
    cutoff_list = get_cutoff_list()
    for cutoff in cutoff_list:
        for method in method_list:
        #for method in ['search_local']:
            row = compare_mapping_algorithms(method_a=method_ref, method_b=method, sim_cutoff_fraction=cutoff)
            row = [cutoff, method_ref, method[0:3]+method[-3:]] + row
            print '\t'.join(map(str,row))

