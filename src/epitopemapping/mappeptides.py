#! /usr/bin/python

import sys
import time
import os


"""
Pseudocode programming process:

GOAL: To accurately map peptides onto proteins.
      If needed, takes into account homologous relationship among source and target proteins.

    (1) map_peptides: Given a list of (A) peptides, (B) target proteins and (C) similarity threshold, 
        retrieve all peptide:protein mappings.
    (2) filter_results: If needed, filter the result to those that meet homologous relationship requirement.
    (3) format_results: Return results.
    
#== Questions:
    Q: Should partial peptide mapping be considered? 
    
#== A list of data validations:
    (1) Warning: Do all source antigens contain non-empty sequences?
    (2) Error: Do all peptides use valid amino acids?
    (3) Stats: How many peptides were mapped? Unmapped? And which peptides?
    (4) Stats: Which proteins had hits? Which didn't?
    (5) Making sure full length of peptides were mapped.  # see blast results.

#== Data Structures
    (1) Results: A text file containing [id_peptide] [id_protein] [similarities] etc.

#== Different techniques to map peptides:
    INPUT: (1) A fasta file containing a list of peptides.
           (2) A fasta file containing a list of source antigens.
           (3) A fasta file containing a list of target antigens.
    method_brute    # Goes over one by one overlapping windows of peptides. Used to create benchmark set.
    method_blast
    method_muscle
    
    
#== Tests ==
    (1) With respect to method_brute,how many mappings are retrieved by methodA?
"""




""" 
Too many input variables. Take out some variables to at most 7 total.
    Ideally:
    INPUT:
        fname_peptide_list_fasta
        fname_source_antigens_blastdb
        fname_target_antigens_blastdb
        d_id_to_peptide
        antigen_list  # Target antigen list; Not sure why this is being used.
        
    == pseudo programming method ==
    1. map peptides to antigens.
    2. If true, keep only best match for each peptide.
    3. If true, keep only those mapping records that involve homologous source and target antigens.
"""

from algorithm.factory_filter_mapping import filter_mapping_unique_keep_ties, filter_mapping_by_homology_protein

report_time = lambda ti,tj, description:  sys.stdout.write(' %8.2f seconds ' %(tj-ti,)+'\t'+description+'\n')

def map_peptides(fname_peptide_list, 
                 fname_antigens_source, 
                 fname_antigens_target, 
                 d_peptides_info, 
                 sim_cutoff_fraction=0.80, 
                 enforce_antigen_homology=False, 
                 filter_mapping=False,
                 method='blast_a',
                 debug=False):
    """
    This version bypasses object creation; Simplify the process of mapping. Separate assigning weights. Make it work with other programs.
        Design001: (1) Each algorithm returns mapping results in same format. (2) hobj used to filter results. (3) Group, reformat results.
        Design002: (1) Each algorithm uses its own mapping results. (2) Same object is used to filter results using hobj.   
    """
    #== 1. Map peptides onto proteins:
    (d_mapping, mappingAlgo) = map_peptides_selector(fname_peptide_list, 
                                                     fname_antigens_target, 
                                                     method=method, 
                                                     sim_cutoff_fraction=sim_cutoff_fraction)

    #== 1. Filter d_mapping to address multiple matches:
    if filter_mapping==True:
        d_mapping = filter_mapping_unique_keep_ties(d_mapping)
    
    #== 2. Filter based on antigen-homology.
    # Calculate homologous relationship between source and target antigens; 
    # Then filtering mapping results based on this homology data.
    ti = time.time()
    hobj = None # Get homology data between source and target antigens
    if enforce_antigen_homology == True:
        (d_mapping, hobj) = filter_mapping_by_homology_protein(d_mapping, d_peptides_info, fname_antigens_source, fname_antigens_target)
    tj = time.time(); td = tj-ti; 
    report_time(ti,tj, '== Map peptides: filter_mapping_by_homology_protein ==')
    
    
    return d_mapping


def map_peptides_selector(fname_peptide_list,  
                          fname_antigens_blastdb_target, 
                          method='blast_a', 
                          sim_cutoff_fraction=1.0):
    """
    Given a peptide:protein mapping algorithm, map a list of peptides onto antigens_target.
    """
    from algorithm.factory import get_mapping_algorithm
    mappingAlgo = get_mapping_algorithm(method=method)
    d_mapping = mappingAlgo.search(fname_peptide_list, fname_antigens_blastdb_target, sim_cutoff_fraction=sim_cutoff_fraction)
    return d_mapping, mappingAlgo
    

def annotate_mapping(d_mapping, antigens_list_target, d_id_to_peptide, debug=False):
    """
    GOAL:
       For each antigen, indicate its mapping status and provide description. 
       (1) Whether there is a peptide mapped to the antigen. 
       (2) Protein description, 
       (3) How many peptides mapped, etc.
       
       d_mapping # Contains mapping records.
       antigens_list_target # A list of all antigens.
       d_id_to_peptide  # For a given peptide id, get its peptide sequence.
       
       How about:
       INPUT:  d_mapping
       OUTPUT: d_annotation  # d[id_antigen] = {description:###, num_peptides:###, mapped:[True,False]}
    """
    from algorithm.blast import parse_sequence_record
    
    d_annotation = {}  # Also include information for unmapped_orfs
    
    for (i, record) in enumerate(antigens_list_target): # For each antigen, return a list of peptides that are found.
        (id_antigen, description) = parse_sequence_record(record)

        match_list = None
        if d_mapping.has_key(id_antigen):
            match_list_temp = d_mapping[id_antigen]
            #match_list = [(sim, dic_id_to_peptide[id_peptide]['linear_peptide_seq'], pos_i,pos_j) for (sim,id_peptide,pos_i,pos_j) in match_list_temp]
            #match_list = [(sim, dic_id_to_peptide[id_peptide]['linear_peptide_seq'], pos_i-1, pos_j-1) for (sim,id_peptide,pos_i,pos_j) in match_list_temp]
            match_list = [(sim, id_peptide, d_id_to_peptide[id_peptide], pos_i-1, pos_j-1) for (sim, id_peptide, pos_i, pos_j) in match_list_temp]
        if match_list != None:
            d_annotation[id_antigen] = {'id_antigen':id_antigen, 'description':description, 'num_peptides':len(match_list), 'mapped':True}
            if debug==True: print '\tmap_peptides', '\t', i, '\t', id_antigen, '\t', record.description
        else:
            d_annotation[id_antigen] = {'id_antigen':id_antigen, 'description':description, 'num_peptides':0, 'mapped':False}
    return d_annotation

def get_d_epitope_id_to_antigen_id(d_mapping):
    """
    Returns a dictionary: d[peptide_id] = a list of (similarity, target_antigen_id)
    """
    d_epitope_id_to_antigen_id = {}

    id_antigen_list = d_mapping.keys(); id_antigen_list.sort()
    for id_antigen in id_antigen_list:
        print 'debug.get_d_epitope_id_to_antigen_id', d_mapping[id_antigen].keys()
        match_list = d_mapping[id_antigen]
        if match_list != None:
            for (similarity, id_peptide, pos_i, pos_j) in match_list:
                if d_epitope_id_to_antigen_id.has_key(id_peptide)==True:
                    temp = d_epitope_id_to_antigen_id[id_peptide]
                    temp.append([similarity, id_antigen])
                    d_epitope_id_to_antigen_id[id_peptide] = temp
                else:
                    temp = [similarity, id_antigen]
                    d_epitope_id_to_antigen_id[id_peptide] = temp
    return d_epitope_id_to_antigen_id


def write_d_epitope_id_to_antigen_id(d_id_to_peptide, d_mapping, dir_prefix='./'):
    """
    GOAL: For each id_peptide, return 
    """
    from util_blastpeptides import write_file, get_line_str
    
    header = ['id_peptide', 'peptide', 'is_mapped']
    content = [header]
    
    id_epitope_list_all = d_id_to_peptide.keys()
    id_epitope_list_all = sorted(id_epitope_list_all, key = lambda id_peptide: int(id_peptide))
    
    dic_epitope_id_to_antigen_id = get_d_epitope_id_to_antigen_id(d_mapping)

    for id_peptide in id_epitope_list_all:
        line = [id_peptide, d_id_to_peptide[id_peptide], dic_epitope_id_to_antigen_id.has_key(id_peptide)]
        content.append(line)
        #print 'check_mapping', '\t', id_peptide, '\t', dic_epitope_id_to_antigen_id.has_key(id_peptide)

    write_file([get_line_str(row,delimiter='\t') for row in content], os.path.join(dir_prefix, 'peptide_list.stats.txt'))
