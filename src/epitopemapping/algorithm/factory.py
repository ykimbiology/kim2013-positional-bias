#! /usr/bin/python

import os


def get_mapping_algorithm(method='bruteforce'):
    from blast      import BlastPeptidesA, BlastPeptidesB
    from bruteforce import BruteForce    
    from usearch import Usearch
    #from experimental.mapping_with_msa import MappingWithMSA
    #from smithwaterman import SmithWaterman
    #from muscle import *
    
    method_list_usearch = ['ublast',  'search_global_fulldp', 'search_local', 'search_local_cmplt', 'usearch_local']
    
    mappingAlgo = None
    if   method == 'bruteforce':        mappingAlgo = BruteForce()
    #elif method == 'bruteforce_cython': mappingAlgo = BruteForceCython()    
    elif method == 'blast_a':           mappingAlgo = BlastPeptidesA()
        
    #== This version used peptide search specific blast command. Returns more hits.
    elif method == 'blast_b':           mappingAlgo = BlastPeptidesB()

    elif method in method_list_usearch: mappingAlgo = Usearch(method=method)
    
    
    #== Experimental ==
    elif method == 'muscle': pass
    
    #== Experimental: This option is for HCV only!!! Uses Multiple sequence alignments of different strains of HCV.
    #         If this is chosen, user-provided peptides will be mapped to **HCV MSA**!!.
    #elif method == 'hcv_msa': mappingAlgo = MappingWithMSA()
    
    #elif method == 'smithwaterman': mappingAlgo = SmithWaterman()
    
    return mappingAlgo


def test_get_mapping_algorithm():
    '''
    d_mapping[id_protein] = a list of 'row'; row = [sim, id_peptide, peptide, pos_i, pos_j]
    '''
    dir_prefix = '/home/yohan/workspace/immunome_browser/src/immunomebrowser.epitopemapping/immunomebrowser/epitopemapping/test_data'
    fname_peptides = os.path.join(dir_prefix, 'peptide_list_small.fasta')
    fname_proteins = os.path.join(dir_prefix, 'protein_list_small.fasta')
    method_name = 'hcv_msa' # bruteforce, hcv_msa
    mappingAlgo = get_mapping_algorithm(method=method_name)
    d_mapping = mappingAlgo.search(fname_peptides, fname_proteins, sim_cutoff_fraction=1.0)
    print d_mapping
    
if __name__ == '__main__':
    test_get_mapping_algorithm()