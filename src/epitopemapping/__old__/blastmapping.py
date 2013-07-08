#! /usr/bin/python

import cPickle
from subprocess import Popen, PIPE

'''Re-name [blastmapping] to [get_homologues.py, to reflect its role.'''

class HomologousObject(object):
    '''Used to capture homologous relationships among proteins.
       Later this information will be used to filter peptide:protein mapping to require that
       only those peptide:protein mappings for homologous source:target antigens will be considered.
       This is because it is possible for a peptide to map to relevant places.'''
    def __init__(self, fname_source_antigen_db, fname_target_db, cmd_blastall='./blastall', debug=False):
        '''Do blast. Store results.'''
        if debug==True: print 'DEBUG running HomologousObject'
        self.cmd_blastall = cmd_blastall
        self.evalue_cutoff = str(0.001)
        self.dic_blast_results = self.run_blast(fname_source_antigen_db, fname_target_db, debug=debug)  # [id_source_antigen] = a list of hits (i.e. id_target_antiegn)

    def run_blast(self,fname_source_antigen_db, fname_target_db, debug=False):

        #self.cmd = [self.cmd_blastall, '-F', 'F', '-p','blastp', '-i', fname_source_antigen_db, '-d', fname_target_db, '-m', '8', '-e', self.evalue_cutoff]
        self.cmd = [self.cmd_blastall, '-F', 'T', '-p','blastp', '-i', fname_source_antigen_db, '-d', fname_target_db, '-m', '8', '-e', self.evalue_cutoff]
        if debug==True: print self.cmd
        #print 'debug cmd ', self.cmd
        output = Popen(self.cmd, stdout=PIPE).communicate()[0]
        
        #== DEBUG:
        #print output
        #for row in output.split('\n'):
        #    print row
        dic_blast_results = self.parse(output, debug=debug)
        return dic_blast_results

    def parse(self, output, debug=False):
        '''For each antigen, list all matching peptides meeting sim threshold.'''
        dic_blast_results = {} # gi --> match_list
        lines=output.split('\n')
        lines = [line for line in lines if len(line.strip()) > 0]
        for line in lines:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            (id_query, id_target, percent_identity, length, t,t, query_i, query_j, target_i,target_j, evalue, bitscore) = line.strip().split()
#            id_query  = id_query.replace('id|','').replace('gi|','')
#            id_target = id_target.replace('id|','').replace('gi|','')
            id_query  = id_query.split('|')[1].strip() # Assumes that id is present in thsi from: gi|####|etc.
            id_target = id_target.split('|')[1].strip()

            #row = (float(percent_identity), id_query,(query_i,query_j),(target_i,target_j))

            '''Update target_i such that it correspond to when query_i=1; for example if query_i=5 instead of 1;'''
            #target_i_int = int(target_i) - (int(query_i) - 1)
           # if target_i_int < 1: target_i_int = 1  # 1 is the first residue.
            #row = (float(percent_identity), id_query, int(target_i), int(target_j))
            row                         = (float(percent_identity), id_query, int(target_i), int(target_j))
            #dic_blast_results[id_query] = dic_blast_results.get(id_target, []) + [id_target]
            dic_blast_results[id_query] = dic_blast_results.get(id_query, []) + [id_target]   # I think this is correct.
            dic_blast_results[id_query] = list(set(dic_blast_results[id_query]))  # Remove redundant ones.
        return dic_blast_results

    def is_homologous_test(self, id_source_antigen, id_target_antigen):
        ''''''
        return True


    #def is_homologous_test(self, id_source_antigen, id_target_antigen):
    def is_homologous(self, id_source_antigen, id_target_antigen):
        '''The two antigens are homologous if their ids are in the hitlist of the former.'''
        #hit_list = self.get_hit_list(id_source_antigen)
        hit_list = []
        if self.dic_blast_results.has_key(id_source_antigen): 
            hit_list = self.dic_blast_results[id_source_antigen]
        return id_target_antigen in hit_list

    def write(self, fname):
        cPickle.dump(self.dic_blast_results, open(fname,'w'))


def test_get_homologues():
    fname_source_antigen_db = '/home/life/workspace/immunome_browser/staticfiles/output_20110926_check_mapping_homology/antigens_source/organism_id_1773.fasta'
    fname_target_db         = '/home/life/workspace/immunome_browser/staticfiles/output_20110926_check_mapping_homology/antigens_target/sequences.1773.fasta'
    cmd_blastall            = '/home/yohan/resource_biology/blast/blast-2.2.23/bin/blastall'
    hobj = HomologousObject(fname_source_antigen_db, fname_target_db, cmd_blastall=cmd_blastall, debug=False)
    hobj.write('hobj.cpickle')
    
    id_source_antigen = '56377'
    id_target_antigen = '57116801'
    
    print '**********'
    result = hobj.dic_blast_results[id_source_antigen]
    for row in result:
        print row
    print '**********'
    
    
    print [id_source_antigen, id_target_antigen, hobj.is_homologous(id_source_antigen, id_target_antigen)]
    assert False
        
if __name__ == '__main__':
    test_get_homologues()
    