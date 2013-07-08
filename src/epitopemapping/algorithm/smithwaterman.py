
'''A wrapper around the smith-waterman algorithm ssearch36.

Binary downloaded from:
http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml

NOTES: 

'''
from subprocess import Popen, PIPE

from main.settings import CMD_SSEARCH

class SmithWaterman(object):
    '''
    TODO Address problem of ssearch not taking longer arguments?
        I want to match full length of the query sequences, rather than subsequence.
    '''
    def __init__(self, debug=False):
        '''Map peptides onto antigens.'''
        self.cmd_ssearch = CMD_SSEARCH
        self.sim_cutoff_percent = 100.0
        
        #== Will hold raw mapping results
        self.output = ''

    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0, debug=False):
        '''public:'''
        self.sim_cutoff_percent = sim_cutoff_fraction*100.0
        
        # Example command: ssearch36 -p peptide_list.fasta sequences.11103.fasta -m 8 -f 100 > output2.txt
        self.cmd = [self.cmd_ssearch, '-p', fname_peptides, fname_proteins, '-m', '8', '-f', '100']
        print ' '.join(self.cmd)
        
        if debug == True: print self.cmd
        self.output = Popen(self.cmd, stdout=PIPE).communicate()[0]
        
        if debug == True:
            for row in self.output.split('\n'): print row
                
        d_mapping = self.parse(self.output, debug=debug)
        return d_mapping

    def get_ncbi_gi(self, id_target):
        #== private:
        ncbi_gi = id_target.split('|')[1].strip()  # Should be a gi.
        assert len(ncbi_gi.strip()) > 0
        return ncbi_gi
        
    def parse(self, output, debug=False):
        '''private: For each antigen, list all matching peptides meeting sim threshold.
        Here, id_target refers to a string version of ncbi_gi number only.'''
        d_mapping = {} # id_target_antigen --> match_list
        lines = output.split('\n')
        
        #== Remove blank lines:
        lines = [line for line in lines if len(line.strip()) > 0]
        
        #== Remove lines with '#' charactr in front.
        lines = [line for line in lines if line[0]!='#']
        
        #== Remove those mappings that do not meet similarity threshold:
        lines_filtered = [line for line in lines if float(line.strip().split()[2]) >= self.sim_cutoff_percent]
        
        for line in lines_filtered:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            
            (id_query, id_target, percent_identity, length, t,t, query_i, query_j, target_i,target_j, evalue, bitscore) = line.strip().split()
            id_query = id_query.replace('id|','')
            #id_target = id_target.replace('gi|','')
            id_target = self.get_ncbi_gi(id_target)

           #I am not sure whether this updating is appropriate:
            '''Update target_i such that it correspond to when query_i=1; for example if query_i=5 instead of 1;'''
            target_i_int = int(target_i) - (int(query_i) - 1)
            if target_i_int < 1: target_i_int = 1  # 1 is the first residue.
            #row = (float(percent_identity), id_query, int(target_i), int(target_j))
            row = (float(percent_identity), id_query, target_i_int, int(target_j))
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


sw = SmithWaterman(debug=True)
