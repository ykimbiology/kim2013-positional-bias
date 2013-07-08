#! /usr/bin/python

'''You can use nosetesting script to run this test.'''

import cPickle
from subprocess import Popen, PIPE

from blastmapping.blastmapping import *

def test_is_homologous():
    '''Test this on the self against self'''
    cmd_blastall = '/home/life/resource_biology/blast/blast-2.2.23/bin/blastall'
    #fname_source_antigen_db = '/home/life/workspace/immunome_browser/data_test/sequences.10245.fasta'
    #fname_source_antigen_db = '/home/life/workspace/immunome_browser/experiments/answer_jon_yewdell/results/organism_5833/query_tcell/data_antigens_target/sequences.5833.fasta'
    #fname_target_antigen_db = fname_source_antigen_db

    # Should be true.
    id_source_antigen = '66275994'
    id_target_antigen = '66275994'

    hobj            = HomologousObject(fname_source_antigen_db, fname_target_antigen_db,cmd_blastall=cmd_blastall,debug=True)
    is_homologous_a = hobj.is_homologous(id_source_antigen, id_target_antigen)
    print 'is_homologous', '\t', '\t'.join(map(str,[is_homologous_a, id_source_antigen, id_target_antigen]))

     # Should be False
    id_source_antigen = '66275995'
    id_target_antigen = '66275994'
    is_homologous_b   = hobj.is_homologous(id_source_antigen, id_target_antigen)
    print 'is_homologous', '\t', '\t'.join(map(str,[is_homologous_b, id_source_antigen, id_target_antigen]))


    #hobj.write_data('dic_blast_results.cpickle')

    #if (is_homologous_a!=True) and (is_homologous_b!=False): assert False
    assert False



if __name__ == '__main__':
    test_is_homologous()