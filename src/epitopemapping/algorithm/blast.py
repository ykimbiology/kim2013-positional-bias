
import cPickle
from subprocess import *
import os
import sys

'''This should be a separate and independent package: blastpeptides
   Should be in a state where one can distribute it.'''


from main.settings import CMD_BLAST_BLASTALL
from epitopemapping.util_blastpeptides import write_file, get_line_str


'''

GOAL: Map peptides onto proteins.

INPUT: (1) peptide_list.fasta, (2) blastdb of target_antigens
OUTPUT: (1) At a given cutoff, list of [epitope_id] [gi] [pos_query] [pos_target]
'''


def parse_sequence_record(record):
    '''Utility function to parse sequence record, a biopython object.'''
    id_antigen  = record.id.split('|')[1].strip()  # Assumes the format: gi|###|XXXXXX
    #description = record.description.split('|')[4]
    description = '|'.join(record.description.split('|')[2:])
    #description = record.description
    seq_antigen = str(record.seq)   
    return (id_antigen, description) 

def get_sample_d_mapping_antigen():
    '''For testing group_by_id_peptide'''
    d_mapping_antigen = {}
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(1.0, 'a', 1, 2)]
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(0.9, 'a', 20, 21)]
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(1.0, 'b', 11, 12)]
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(1.0, 'b', 21, 22)]
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(0.5, 'b', 41, 42)]
    d_mapping_antigen['1'] = d_mapping_antigen.get('1',[]) + [(0.9, 'b', 31, 32)]
    return d_mapping_antigen


class BlastPeptidesA(object):
    def __init__(self, debug=False):
        '''Map peptides onto antigens.'''
        self.cmd_blastall = CMD_BLAST_BLASTALL
        self.sim_cutoff_percent = 100.0
        
        #== Will hold raw mapping results
        self.output = ''

    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0, debug=False):
        '''public:'''
        self.sim_cutoff_percent = sim_cutoff_fraction*100.0
        
        self.cmd = [self.cmd_blastall, '-F', 'F', '-p', 'blastp', '-m', '8', '-i', fname_peptides, '-d', fname_proteins]
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
        
        #== Remove those mappings that do not meet similarity threshold:
        lines_filtered = [line for line in lines if float(line.strip().split()[2]) >= self.sim_cutoff_percent]
        
        for line in lines_filtered:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            #print 'DEBUG blast', line
            
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


class BlastPeptidesC(object):
    '''This version is for BLAST 2.2.26+ '''
    def __init__(self, debug=False):
        '''Map peptides onto antigens.'''
        self.cmd_blastall = 'blastp'
        self.sim_cutoff_percent = 100.0
        
        #== Will hold raw mapping results
        self.output = ''

    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0, debug=False):
        '''public:'''
        self.sim_cutoff_percent = sim_cutoff_fraction*100.0
        
        # nohup blastp -task blastp-short -evalue 20000.0 -num_threads 5-outfmt 6 -query peptide_list.hla-b-0702.eluted.fasta -db /home/yohan/resource_biology/blast_db/nr/nr -out blast_output.nr.hla-b-0702.eluted > log.txt
        #self.cmd = [self.cmd_blastall, '-task', 'blastp-short',  '-evalue', '10.0', '-num_threads', '5', '-outfmt', '6', '-query', fname_peptides, '-db', fname_proteins]
        #self.cmd = [self.cmd_blastall, '-task', 'blastp-short',  '-evalue', '100.0', '-num_threads', '5', '-outfmt', '6', '-query', fname_peptides, '-db', fname_proteins]
        #self.cmd = [self.cmd_blastall, '-task', 'blastp-short',  '-evalue', '1000.0', '-num_threads', '5', '-outfmt', '6', '-query', fname_peptides, '-db', fname_proteins]
        #self.cmd = [self.cmd_blastall, '-task', 'blastp-short',  '-evalue', '10000.0', '-num_threads', '5', '-outfmt', '6', '-query', fname_peptides, '-db', fname_proteins]
        self.cmd = [self.cmd_blastall, '-task', 'blastp-short',  '-evalue', '20000.0', '-num_threads', '5', '-outfmt', '6', '-query', fname_peptides, '-db', fname_proteins]
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
        
        #== Remove those mappings that do not meet similarity threshold:
        lines_filtered = [line for line in lines if float(line.strip().split()[2]) >= self.sim_cutoff_percent]
        
        for line in lines_filtered:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            #print 'DEBUG blast', line
            
            (id_query, id_target, percent_identity, length, t,t, query_i, query_j, target_i,target_j, evalue, bitscore) = line.strip().split()
            id_query = id_query.replace('id|','')
            #id_target = id_target.replace('gi|','')
            id_target = self.get_ncbi_gi(id_target)

           #I am not sure whether this updating is appropriate:
            #'''Update target_i such that it correspond to when query_i=1; for example if query_i=5 instead of 1;'''
            #target_i_int = int(target_i) - (int(query_i) - 1)
            #if target_i_int < 1: target_i_int = 1  # 1 is the first residue.
            row = (float(percent_identity), id_query, int(target_i), int(target_j))
            #row = (float(percent_identity), id_query, target_i_int, int(target_j))
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
        
def get_peptide_region_length(line):
    index_q_start = 6     # query start
    index_q_end   = 7
    row = line.strip().split()
    q_start = row[index_q_start].strip()
    q_end   = row[index_q_end].strip()
    peptide_region_length = int(q_end) - int(q_start) + 1
    return peptide_region_length
    
def filter_blast_output_table(blast_output, percent_sim_cutoff=80.0):
    index_percent_id = 2

    lines = blast_output.split('\n')
    lines = [line for line in lines if len(line.strip()) > 0]
    
    '''#== Filter Rule #1:'''
    lines_filtered_a = [line for line in lines if float(line.strip().split()[index_percent_id]) >= percent_sim_cutoff]
    
    '''#== Filter Rule #2:'''
    min_peptide_region_length = 7
    lines_filtered_b = []
    for line in lines_filtered_a:
        pregion_length = get_peptide_region_length(line)
        if pregion_length >= min_peptide_region_length:
            lines_filtered_b.append(line)
            
    '''#== Filter Rule #3: For each peptide:protein mapping, (a) assume that 100% match was observed and (b) remove if fraction of peptide_length is less than the similarity threshold:'''
    
    
    return lines_filtered_b
        

class BlastPeptidesB(object):
    '''For [hcv, bcell], blast_b is also not picking up some peptides at 30% sequence identity.
       This version uses peptide-specific blast command. Note the difference in the command used.
       This should eventually replace [BlastPeptide], because it will most likely pick up more valid peptide:protein mappings.
       Blast command taken from: http://www.ncbi.nlm.nih.gov/blast/Why.shtml
       NOTE: BLAST residue positions use 1-index system. That is, first residue gets an index of 1 rather than 0.'''
    def __init__(self, debug=False):
        '''Maps peptides onto antigens.'''
        self.cmd_blastall = CMD_BLAST_BLASTALL
        self.percent_sim_cutoff = 1.0*100
        
        #== Will hold blast output
        self.output = ''
        self.output_f = '' # filtered blast output table.
        
    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0, debug=False):
        self.percent_sim_cutoff = sim_cutoff_fraction*100
        
        self.cmd = [self.cmd_blastall, '-p', 'blastp', '-M', 'PAM30', '-W', '2', '-F', 'F', '-C', 'F', '-m', '8', '-e', '20000', '-i', fname_peptides, '-d', fname_proteins]
        if debug==True: print self.cmd
        self.output = Popen(self.cmd, stdout=PIPE).communicate()[0]
        if debug == True:
            for row in self.output.split('\n'):
                print row
        d_blast_result = self.parse(self.output, debug=debug)
        return d_blast_result

    def parse(self, output, debug=False):
        '''For each antigen, list all matching peptides meeting sim threshold.'''
        lines_filtered = filter_blast_output_table(output, percent_sim_cutoff=self.percent_sim_cutoff)
        
        self.output_f = '\n'.join(lines_filtered)
        
        d_blast_result = {} # gi --> match_list
        for line in lines_filtered:
            #'''id|1446    gi|66275888    100.00    40    0    0    1    40    31    70    3e-20    85.1'''
            if debug==True: print line
            (id_query, id_target, percent_identity, length, t,t, query_i, query_j, target_i,target_j, evalue, bitscore) = line.strip().split()
            id_query = id_query.replace('id|','')
            id_target = id_target.replace('gi|','')
            #row = (float(percent_identity), id_query,(query_i,query_j),(target_i,target_j))

            '''Update target_i such that it correspond to when query_i=1; for example if query_i=5 instead of 1;
            It may be that there are gaps/deletions, thereby resulting in mismatch in two lengths.
            ToDo: Address this problem thoroughly.'''
            target_i_int = int(target_i) - (int(query_i) - 1)
            if target_i_int < 1: target_i_int = 1  # 1 is the first residue.
            #row = (float(percent_identity), id_query, int(target_i), int(target_j))
            row = (float(percent_identity), id_query, target_i_int, int(target_j))
            d_blast_result[id_target] = d_blast_result.get(id_target, []) + [row]
        return d_blast_result
    
    def write(self, fname_output):
        '''Write out raw blast output for viewing later.'''
        f=open(fname_output,'w')
        f.write(self.output)
        # DEBUG:
        print self.output
        f.close()

        #== Writing out filtered blast output table:
        f=open(fname_output.replace('.txt', '.filtered.txt'),'w')
        f.write(self.output_f)
        # DEBUG:
        print self.output_f
        f.close()


class HomologousObject(object):
    '''Used to capture homologous relationships among proteins.
       Later this information will be used to filter peptide:protein mapping to require that
       only those peptide:protein mappings for homologous source:target antigens will be considered.
       This is because it is possible for a peptide to map to irrelevant places.'''
    def __init__(self, fname_source_antigen_db, fname_target_db, debug=False):
        '''Do blast. Store results.'''
        if debug==True: print 'DEBUG running HomologousObject'
        self.cmd_blastall = CMD_BLAST_BLASTALL
        self.evalue_cutoff = str(0.001)
        self.dic_blast_results = self.run_blast(fname_source_antigen_db, fname_target_db, debug=debug)  # [id_source_antigen] = a list of hits (i.e. id_target_antiegn)

    def run_blast(self, fname_source_antigen_db, fname_target_db, debug=False):
        '''Parameter description:
           -F Filter query sequence (DUST with blastn, SEG with others) [String]; default True '''
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

    def check_homologous(self, id_source_antigen, id_target_antigen):
        '''The two antigens are homologous if their ids are in the hitlist of the former.'''
        #hit_list = self.get_hit_list(id_source_antigen)
        hit_list = []
        if self.dic_blast_results.has_key(id_source_antigen): 
            hit_list = self.dic_blast_results[id_source_antigen]
        return id_target_antigen in hit_list
    
    def check_homologous_list(hobj, id_source_antigen_list, id_target_antigen):
        '''Loops over the list of id_source and returns true if any combination is homologous.'''
        is_homologous = False
        for id_source in id_source_antigen_list:
            if hobj.check_homologous(id_source, id_target_antigen)==True:
                is_homologous = True
                break
        return is_homologous

    def write(self, fname):
        cPickle.dump(self.dic_blast_results, open(fname,'w'))


def test_get_homologues():
    '''ToDo: Need to update the paths to input files.
    '''
    fname_antigen_source    = '/home/yohan/workspace/immunome_browser/staticfiles/output_20110926_check_mapping_homology/antigens_source/organism_id_1773.fasta'
    fname_antigen_target    = '/home/yohan/workspace/immunome_browser/staticfiles/output_20110926_check_mapping_homology/antigens_target/sequences.1773.fasta'
    cmd_blastall            = '/home/yohan/resource_biology/blast/blast-2.2.23/bin/blastall'

    hobj = HomologousObject(fname_antigen_source, fname_antigen_target, cmd_blastall=cmd_blastall, debug=False)
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
    pass
    #test()