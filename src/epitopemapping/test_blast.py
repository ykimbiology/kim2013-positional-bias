import os
import cPickle
import sys


from main.settings import CMD_BLAST_BLASTALL
from main.util_antigens import *

cmd_blastall = CMD_BLAST_BLASTALL

# Local:
#from blast import map_peptides_blast

#def read_data():
#    '''Preparing data for testing:'''
#    dir_prefix = '/home/life/workspace/immunome_browser/staticfiles/output_20110926_check_mapping_homology'
#    #dir_prefix = '/home/life/workspace/immunome_browser/data_test/epitope_antigen_mapping/'
#    fname_dic_output              = os.path.join(dir_prefix, 'plotting', 'data_mapping.cpickle')
#    fname_peptide_list_fasta      = os.path.join(dir_prefix, 'peptides_info', 'peptide_list.fasta')
#    fname_source_antigens_blastdb = os.path.join(dir_prefix, 'antigens_source', 'organism_id_1773.fasta')
#    fname_target_antigens_blastdb = os.path.join(dir_prefix, 'antigens_target', 'sequences.1773.fasta')
#
#    data_temp = cPickle.load(open(fname_dic_output, 'r'))
#
#    d_output = {}
#    #fraction_sim_cutoff = 0.40
#    fraction_sim_cutoff = 0.80
#    d_output['fraction_sim_cutoff']           = fraction_sim_cutoff
#    d_output['fname_peptide_list_fasta']      = fname_peptide_list_fasta
#    d_output['fname_source_antigens_blastdb'] = fname_source_antigens_blastdb
#    d_output['fname_target_antigens_blastdb'] = fname_target_antigens_blastdb
#    #cPickle.dump([dic_peptide_info, dic_id_to_peptide, dic_antigens, dic_peptide_weight, dic_content_stats], fout_plotting)
#    d_output['dic_peptide_info']   = data_temp[0]
#    d_output['dic_id_to_peptide']  = data_temp[1]
#    d_output['dic_antigens']       = data_temp[2]
#    d_output['dic_peptide_weight'] = data_temp[3]
#    d_output['dic_content_stats']  = data_temp[4]
#    return d_output

    #dic_content_stats = map_peptides_blast(fname_peptide_list_fasta, fname_target_antigens_blastdb, dic_id_to_peptide, antigen_list, fraction_sim_cutoff=fraction_sim_cutoff, cmd_blastall = cmd_blastall, debug=False)
    #dic_content_stats = d['dic_content_stats']

    #id_antigen_list = dic_content_stats.keys()
    #for id_antigen in id_antigen_list:
    #    print id_antigen, dic_content_stats[id_antigen]['mapped'], dic_content_stats[id_antigen]['num_peptides']

    #assert False


#
#def test():
#    dir_prefix = '/home/life/workspace/immunome_browser/experiments/20100812_optimize_mapping/data'
#    name_peptide_list = os.path.join(dir_prefix,'peptide_list.fasta')
#    name_db = os.path.join(dir_prefix,'sequences.10245.fasta')
#    dir_blast = '/home/life/resource_biology/blast/blast-2.2.23/bin'
#
#    blast = BlastPeptide(dir_blast=dir_blast)
#    dic = blast.search(name_peptide_list, name_db)
#    print 'len ', len(dic)
#    key_list = dic.keys()
#    for key in key_list:
#        print key, '\t', dic[key]
#
#def test_map_peptides():
##    dir_prefix = '/home/life/workspace/immunome_browser/experiments/20100812_optimize_mapping/data'
##    fname_peptide_list_fasta = os.path.join(dir_prefix,'peptide_list.fasta')
##    fname_target_antigens_blastdb = os.path.join(dir_prefix,'sequences.10245.fasta')
#
#
#    fname_peptide_list_fasta = os.path.join('/home/life/workspace/immunome_browser/experiments/immunome_browser_core/www/output_1773/data_peptides','peptide_list.fasta')
#    fname_target_antigens_blastdb = os.path.join('/home/life/workspace/immunome_browser/experiments/immunome_browser_core/www/output_1773/data_target_antigens','sequences.1773.fasta')
#
#    dic_homology = None
#    (dic_id_to_peptide, dic_antigens) = read_data()
#    dic_content_stats = map_peptides_blast(fname_peptide_list_fasta, fname_target_antigens_blastdb, dic_id_to_peptide, dic_antigens.values(), dic_homology, fraction_sim_cutoff=0.80, debug=False)
#    id_antigen_list = dic_content_stats.keys()
#    for id_antigen in id_antigen_list:
#        print id_antigen, dic_content_stats[id_antigen]['mapped'], dic_content_stats[id_antigen]['num_peptides']
#
#
#    print '== id_epitope --> [list of id_antigen] =='
#    id_epitope_list_all = dic_id_to_peptide.keys(); id_epitope_list_all.sort()
#
#    dic_epitope_id_to_antigen_id = get_dic_epitope_id_to_antigen_id(dic_content_stats)
#    #id_epitope_list = dic_epitope_id_to_antigen_id.keys()
#    #for id_epitope in id_epitope_list:
#    #    temp = dic_epitope_id_to_antigen_id[id_epitope]
#    #    print id_epitope, '\t', len(temp), '\t', temp
#    for id_epitope in id_epitope_list_all:
#        print 'check_mapping', '\t', id_epitope, '\t', dic_epitope_id_to_antigen_id.has_key(id_epitope)




def test_is_homologous():
    '''Test this on the self against self'''
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
    (dic_content_stats, blast, hobj) = test_map_peptides_blast()

    #d = read_data()
    #print 'len(d)', len(d), d.keys()