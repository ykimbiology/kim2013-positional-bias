#! /usr/bin/python

import os
import sys
sys.path.append('/home/yohan/workspace/immunomebrowser/core/src/immunomebrowser.main/immunomebrowser')


path_prefix = '/home/yohan/data_large/immunomebrowser/experiments/test_data'
fname_peptide_list = os.path.join(path_prefix, 'peptide_list.fasta')
fname_protein_list = os.path.join(path_prefix, 'sequences.10245.fasta')

from usearch import Usearch

def test_usearch():
    #method_list = ['ublast', 'search_global', 'search_global_fulldp', 'search_local', 'usearch_local', 'usearch_global']
    method_list = ['ublast',  'search_global_fulldp', 'search_local', 'usearch_local']
    for method in method_list:
        malgo = Usearch(method=method)
        d_mapping = malgo.search(fname_peptide_list, fname_protein_list, sim_cutoff_fraction=1.0)
        print '=======', method, len(d_mapping), '============'
        #id_antigen_list = d_mapping.keys()
        #for id_antigen in id_antigen_list:
            #print [id_antigen, len(d_mapping[id_antigen])]



test_usearch()