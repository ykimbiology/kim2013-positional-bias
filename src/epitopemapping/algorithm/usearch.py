
"""
A wrapper around the usearch algorithm.

Binary downloaded from:
http://www.drive5.com/usearch

NOTES: 

"""

import os
import random
from subprocess import Popen, PIPE

from main.settings import CMD_USEARCH

def get_method_list():
    method_list = ['ublast',  'search_global_fulldp', 'search_local', 'usearch_local']
    return method_list

class Usearch(object):
    """
    Does global alignment rather than local.
    One can specify either (1) heuristic or (2) full dynamic programming.
    
    == Different combinations of commands ==
    -search_global -id [###]                  #search_global
    -search_global -id [###] -fulldp          #search_global_fulldp
    -search_local -evalue 10.0                #search_local
    -ublast -evalue 10.0                      #ublast
    -usearch_global -id [###]                 #usearch_global
    -usearch_local -evalue 10.0               #usearch_local
    
    """
    def __init__(self, method='ublast', debug=False):
        '''Map peptides onto antigens.'''
        self.cmd_usearch = CMD_USEARCH
        self.sim_cutoff_percent = 100.0
        self.method_type = method
        self.cmd = []  # initializing to empty set of arguments.
        
        #== Will hold raw mapping results
        self.output = ''

    def get_d_cmd(self, fname_peptides, fname_proteins, fname_output, sim_cutoff_fraction):
        # Local alignment options not supported by -ublast
        d_cmd={}
        d_cmd['ublast']               = [self.cmd_usearch, '-ublast', fname_peptides, '-db', fname_proteins, '-evalue', '10.0', '-query_cov', '1.0', '-blast6out', fname_output]
        d_cmd['search_global']        = [self.cmd_usearch, '-search_global', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-query_cov', '1.0', '-blast6out', fname_output]
        d_cmd['search_global_fulldp'] = [self.cmd_usearch, '-search_global', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-query_cov', '1.0', '-fulldp', '-blast6out', fname_output]
        d_cmd['search_local']         = [self.cmd_usearch, '-search_local', fname_peptides, '-db', fname_proteins, '-evalue', '10.0', '-blast6out', fname_output]
        #d_cmd['search_local_cmplt']   = [self.cmd_usearch, '-search_local', fname_peptides, '-db', fname_proteins, '-evalue', '10.0', '-query_cov', '1.0', '-blast6out', fname_output]
        d_cmd['search_local_cmplt']   = [self.cmd_usearch, '-search_local', fname_peptides, '-db', fname_proteins, '-evalue', '10.0', '-query_cov', '1.0', '-lopen', '100.0', '-lext', '1.0', '-blast6out', fname_output]
        d_cmd['usearch_local']        = [self.cmd_usearch, '-usearch_local', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-query_cov', '1.0', '-lopen', '100.0', '-lext', '1.0', '-blast6out', fname_output]
        d_cmd['usearch_global']       = [self.cmd_usearch, '-usearch_global', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-blast6out', fname_output]
        
        # For tb, only 4 hits.
        #d_cmd['search_global']        = [self.cmd_usearch, '-search_global', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-blast6out', fname_output]        
        # For tb, only 4 hits.
        #d_cmd['usearch_global']       = [self.cmd_usearch, '-usearch_global', fname_peptides, '-db', fname_proteins, '-id', sim_cutoff_fraction, '-blast6out', fname_output]
        return d_cmd
    
    def get_fname_output(self):
        '''Randomly generate a file name to contain the search results. '''
        fname_output = '.'.join(['output', str(random.random())[2:], 'txt'])
        return fname_output
    
    def get_path_dir(self, fname):
        row = fname.split('/')
        path_dir = '/'.join(row[0:-1])
        return path_dir
    
    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0, debug=False):
        '''public:'''
        
        self.sim_cutoff_percent = sim_cutoff_fraction*100.0
        path_dir = self.get_path_dir(fname_peptides)
        fname_output = self.get_fname_output()
        fname_output = os.path.join(path_dir, fname_output)
        
        self.d_cmd = self.get_d_cmd(fname_peptides, fname_proteins, fname_output, sim_cutoff_fraction)
        
        #Example: usearch6.0.152_i86linux32  -search_global peptide_list.10k.fasta -db sequences.10245.fasta  -id 0.5 -fulldp -alnout output.txt
        
        self.cmd = self.d_cmd[self.method_type]
        self.cmd = map(str, self.cmd)
        print ' '.join(self.cmd)
        
        if debug == True: print self.cmd
        self.output = Popen(self.cmd, stdout=PIPE).communicate()[0]
        
        if debug == True:
            for row in self.output.split('\n'): print row
                
        d_mapping = self.parse(fname_output, debug=debug)
        # Remove the output file:
        os.remove(fname_output)
        
        return d_mapping

    def get_ncbi_gi(self, id_target):
        #== private:
        ncbi_gi = id_target.split('|')[1].strip()  # ex: gi|###, id|###
        assert len(ncbi_gi.strip()) > 0
        return ncbi_gi
        
    def parse(self, fname_output, debug=False):
        '''private: For each antigen, list all matching peptides meeting sim threshold.
        Here, id_target refers to a string version of ncbi_gi number only.'''
        d_mapping = {} # id_target_antigen --> match_list
        
        lines = open(fname_output,'r').read().split('\n')
        
        #== Remove blank lines:
        lines = [line for line in lines if len(line.strip()) > 0]
        
        #== Remove those mappings that do not meet similarity threshold:
        lines_filtered = [line for line in lines if float(line.strip().split()[2]) >= self.sim_cutoff_percent]
        
        for line in lines_filtered:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            
            (id_query, id_target, percent_identity, length, t,t, query_i, query_j, target_i,target_j, evalue, bitscore) = line.strip().split()
            id_query = id_query.replace('id|','')
            #id_target = id_target.replace('gi|','')
            id_target = self.get_ncbi_gi(id_target)

            '''Positions in 1-index'''
            #row = (float(percent_identity), id_query, int(target_i), int(target_j))
            row = (float(percent_identity), id_query, int(target_i), int(target_j))
            d_mapping[id_target] = d_mapping.get(id_target, []) + [row]
            
            if debug==True: id_target, row
        return d_mapping
    
    def write(self, fname_output):
        '''public: Write out raw blast output for viewing later.'''
        f=open(fname_output,'w')
        f.write(self.output)
        # DEBUG:
        #print self.output
        f.close()


