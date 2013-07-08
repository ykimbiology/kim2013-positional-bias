#! /usr/bin/python


import os
import sys
import cPickle

import Bio.AlignIO

'''
== Algorithm ==
    1. For each polyprotein, generate 'd_mapping'
    2. Combine a list of 'd_mapping'.
    3. Group by 'id_peptide'.
    4. For each peptide, only keep best matches.
    5. Convert it into 'd_mapping_resolved'. # Mapped positions are with respect to 'H77'.
'''

sys.path.append('/home/yohan/workspace/immunome_browser/src/immunomebrowser.epitopemapping/immunomebrowser')
from epitopemapping.algorithm.factory import get_mapping_algorithm


#== Settings =================================================================
DIR_PREFIX = '/home/yohan/workspace/immunome_browser/experiments/hcv/20111129_hcv_mapping_stats/data'
DIR_PREFIX_SEQUENCES = '/home/yohan/workspace/immunome_browser/experiments/hcv/20111129_hcv_mapping_stats/sequences'
FNAME_PEPTIDE_LIST_FASTA = os.path.join(DIR_PREFIX, 'peptide_list.fasta') #== Will be variable.
FNAME_PROTEINS = os.path.join(DIR_PREFIX_SEQUENCES, 'all.faa')            #== Fixed.
FNAME_MSA      = os.path.join(DIR_PREFIX_SEQUENCES, 'msa.aln')            #== Fixed. msa_sample.aln, msa.aln


class PositionMappingWithMSA(object):
    '''
    GOAL: Convert positions of peptides mapped to one protein to those of the reference protein.
    
    pos_mapping = PositionMapingWithMSA()
    pos_mapping.get_position_for_msa(id_protein, pos_i, pos_j)
    pos_ref = pos_mapping.to_protein_reference(id_protein_ref, id_protein, pos_i, pos_j)
    '''
    def __init__(self):
        self.char_gap = '-'   # Character for the gaps:
        (d_msa, d_seq) = self.read_msa(FNAME_MSA)
        self.d_pos_mapping = self.get_d_pos_mapping(d_seq, d_msa)
        
    def parse_id_record(self, id_record):
        ''' 
        INPUT: gi|1234|abcdef
        OUTPUT: 1234, in string format.
        '''
        id_record_parsed = id_record.split('|')[1].strip()
        return id_record_parsed
    
    def test_get_position_msa(self):
        seq = 'ABCDEF'
        seq_msa = 'AB-CDE-F'
        #pos = 3
        pos = 6
        pos_msa = self.get_position_msa(seq, seq_msa, pos)
        print seq
        print seq_msa
        print [pos, pos_msa]
        
    def get_position_msa(self, seq, seq_msa, pos):
        '''Returns position of the sequence in the context of msa.
        '''
        index = pos - 1   # So that first residue is zero.
        count_res = 0   # Counts number of residues, excluding gaps.
        index_msa = None
        for i,res in enumerate(seq_msa):
            if res != self.char_gap: count_res = count_res + 1
            if count_res == pos:
                index_msa = i
                break
        pos_msa = index_msa + 1  # Residue starts at one.
        return pos_msa
            
    def get_d_pos_seq_to_msa(self, seq, seq_msa):
        '''
        INPUT:
            seq = sequence without gaps.
            seq_msa = sequence with gaps, as stored in the multiple seq alignment.
        Returns a dictionary that maps a position from a sequence to that of msa.
        '''
        d_pos_seq_to_msa = {}
        d_pos_msa_to_seq = {}
        for i,res in enumerate(seq):
            pos_seq = i + 1
            pos_msa = self.get_position_msa(seq, seq_msa, pos_seq)
            d_pos_seq_to_msa[pos_seq] = pos_msa
            d_pos_msa_to_seq[pos_msa] = pos_seq
        return (d_pos_seq_to_msa, d_pos_msa_to_seq)
    
    # Looks like in terms of speed, not much improvement over [get_d_pos_seq_to_msa].
    def get_d_pos_seq_to_msa_fast(self, seq, seq_msa):
        '''
        A faster version.
        INPUT:
            seq = sequence without gaps.
            seq_msa = sequence with gaps, as stored in the multiple seq alignment.
        Returns a dictionary that maps a position from a sequence to that of msa.
        '''
        d_pos_seq_to_msa = {}
        d_pos_msa_to_seq = {}
        count_seq = 0 # Number of valid residues counted in the input sequence.
        count_msa = 0 # Number of valid residues counted in the corresponding seq in msa with gaps.
        index_msa = 0
        for index_seq, res in enumerate(seq):
            count_seq = count_seq + 1
            while count_msa != count_seq:
                res_msa = seq_msa[index_msa]
                index_msa = index_msa + 1
                if res_msa != self.char_gap: count_msa = count_msa + 1
            pos_seq = index_seq + 1
            pos_msa = index_msa + 1
            d_pos_seq_to_msa[pos_seq] = pos_msa
            d_pos_msa_to_seq[pos_msa] = pos_seq
        return (d_pos_seq_to_msa, d_pos_msa_to_seq)
    
    def get_d_pos_mapping(self, d_seq, d_msa):
        d_pos_mapping = {}   # d_pos_mapping[id_protein] = (d_pos_seq_to_msa, d_pos_msa_to_seq)
        id_protein_list = d_seq.keys()
        for id_protein in id_protein_list:
            seq = d_seq[id_protein]
            seq_msa = d_msa[id_protein]
            (d_pos_seq_to_msa, d_pos_msa_to_seq) = self.get_d_pos_seq_to_msa(seq, seq_msa)
            d_pos_mapping[id_protein] = (d_pos_seq_to_msa, d_pos_msa_to_seq)
        return d_pos_mapping
    
    def read_msa(self, fname_msa):
        d_msa = {}  # sequences still containing the gap characters.
        d_seq = {}  # sequences after removing gap characters.   
        msa = Bio.AlignIO.read(open(fname_msa), 'fasta')
        for record in msa:
            id_record = self.parse_id_record(record.id)
            seq = str(record.seq)
            d_msa[id_record] = seq
            d_seq[id_record] = seq.replace(self.char_gap, '')
        return d_msa, d_seq

    def to_msa(self, id_protein, pos_i, pos_j):
        '''
        INPUT: (pos_i, pos_j): Residue positions start at one, rather than zero.
        OUTPUT: (pos_i_msa, pos_j_msa)
        '''
        (d_pos_seq_to_msa,d_pos_msa_to_seq) = self.d_pos_mapping[id_protein]
        pos_i_msa = d_pos_seq_to_msa[pos_i]
        pos_j_msa = d_pos_seq_to_msa[pos_j]
        return (pos_i_msa, pos_j_msa)
    
    def to_protein_reference(self, id_protein_ref, id_protein, pos_i, pos_j):
        '''
        Convert from positions in one protein to that of another.
        None returned if deletion in the target protein, according to the msa.
        '''
        
        #== 1. Get positions with respect to msa:
        (pos_i_msa, pos_j_msa) = self.to_msa(id_protein, pos_i, pos_j)
        
        #== 2. Get positions with respect to the protein_ref:
        (d_pos_seq_to_msa, d_pos_msa_to_seq) = self.d_pos_mapping[id_protein_ref]
        pos_i_ref = None
        pos_j_ref = None
        if d_pos_msa_to_seq.has_key(pos_i_msa)==True:
            pos_i_ref = d_pos_msa_to_seq[pos_i_msa]
        if d_pos_msa_to_seq.has_key(pos_j_msa)==True:    
            pos_j_ref = d_pos_msa_to_seq[pos_j_msa]
        return (pos_i_ref, pos_j_ref)




class MappingWithMSA(object):
    '''
    GOAL: Map peptides to antigen using multiple sequence alignment.
    Presumably using MSA will result in more reliable peptide:antigen mapping.
    This should result in more peptides with higher percent sequence identities.
    '''
    def __init__(self, method='bruteforce', debug=False):
    #def __init__(self, assay_type='bcell', mode='generate_data', method='bruteforce', debug=False):
        '''
        mode=[generate_data = run mapping; otherwise read mapping results from file.]
        '''
        self.debug = debug
        self.method = method  # Specify method for mapping peptides to proteins.
        self.output = {}    # Store intermediate data for debugging later.
    
    def search(self, fname_peptide_list_fasta, fname_target_antigens_blastdb, sim_cutoff_fraction=0.80):
        '''
        This version tries to address presence of those peptides that could not be mapped.
        NOTE: fname_target_antigens_blastdb is not used, b/c peptides are mapped to the list of HCV strain polyproteins.
        Uses multiple sequence alignments to map peptides to antigens.
        OUTPUT: 
            d_mapping[id_protein] = a list of 'row', where row = [sim_percent, id_peptide, pos_protein_i, pos_protein_j]
        '''
        
        mappingAlgo = get_mapping_algorithm(method=self.method)
        
        #== Map all peptides to the list of polyproteins.
        d_mapping = self.get_d_mapping(fname_peptide_list_fasta, sim_cutoff_fraction, mappingAlgo)
        #d_mapping = self.get_d_mapping_sample()
     
        #== Map all positions to those of H77; Remove any that overlap with deleted region in H77.
        d_mapping_resolved = self.resolve_d_mapping(d_mapping)
             
     
        match_list_all = self.combine_mapping_list(d_mapping_resolved)
        d_idpeptide = self.group_by_id_peptide(match_list_all)
        
        #== For each peptide, keep only the best hits. There may be duplicates.
        d_idpeptide = self.filter_mapping_unique_keep_ties(d_idpeptide)
        match_list_all_f = d_idpeptide.values()
        match_list_all_after = []
        [match_list_all_after.extend(row) for row in match_list_all_f]
        
        d_mapping_after = self.group_by_id_protein(match_list_all_after)
        
        #== Storing intermediate data:
        self.output['d_mapping']          = d_mapping
        self.output['d_mapping_resolved'] = d_mapping_resolved
        self.output['d_mapping_after']    = d_mapping_after
        
        

        return d_mapping_after
    
    def search_old(self, fname_peptide_list_fasta, fname_target_antigens_blastdb, sim_cutoff_fraction=0.80):
        '''
        NOTE: fname_target_antigens_blastdb is not used, b/c peptides are mapped to the list of HCV strain polyproteins.
        Uses multiple sequence alignments to map peptides to antigens.
        OUTPUT: 
            d_mapping[id_protein] = a list of 'row', where row = [sim_percent, id_peptide, pos_protein_i, pos_protein_j]
        '''
        
        mappingAlgo = get_mapping_algorithm(method=self.method)
        d_mapping = self.get_d_mapping(fname_peptide_list_fasta, sim_cutoff_fraction, mappingAlgo)
        #d_mapping = self.get_d_mapping_sample()
     
        match_list_all = self.combine_mapping_list(d_mapping)
        d_idpeptide = self.group_by_id_peptide(match_list_all)
        
        d_idpeptide = self.filter_mapping_unique_keep_ties(d_idpeptide)
        match_list_all_f = d_idpeptide.values()
        match_list_all_after = []
        [match_list_all_after.extend(row) for row in match_list_all_f]
        
        d_mapping_after = self.group_by_id_protein(match_list_all_after)
        
        #d_mapping_resolved = d_mapping_after   #== Uncomment if you want to skip the 'resolving' step.
        d_mapping_resolved = self.resolve_d_mapping(d_mapping_after)

        return d_mapping_resolved
    
    
    def resolve_d_mapping(self, d_mapping):
        ''' All mappings will be using residue numbering with respect to the reference protein.
            Resolve duplicate matches by processing [d_mapping_after].
            This involves:
            For non-reference protein, translating positions to those of reference protein. 
        '''
        print 'debug: resolve_d_mapping, d_mapping', d_mapping
        pos_mapping = PositionMappingWithMSA()
        
        id_protein_ref = '22129793'  # HCV H77 sequence.
        id_protein_list = d_mapping.keys()
        id_protein_list_f = [id for id in id_protein_list if id != id_protein_ref]
        match_list_all = []
        match_list_all = match_list_all + d_mapping.get(id_protein_ref, [])
        match_list_none = []   #=A list of those matches that mapped to gapped region in ref protein.
        
        d_mapping_resolved = {}
        for id_protein in id_protein_list_f:
            match_list_temp = d_mapping[id_protein]
            match_list_translated = []
            for match in match_list_temp:
                (similarity, id_peptide, pos_i, pos_j) = match
                (pos_i_ref, pos_j_ref) = pos_mapping.to_protein_reference(id_protein_ref, id_protein, pos_i, pos_j)
                match_resolved = (similarity, id_peptide, pos_i_ref, pos_j_ref)
                
                '''Do not add the match if any one position is None, and thus mapped to gapped region in the reference protein.'''
                if (pos_i_ref == None) or (pos_j_ref == None): 
                    match_list_none.append((id_protein, match, match_resolved))
                else: 
                    match_list_all.append(match_resolved)
                
                
        #d_mapping_resolved[id_protein_ref] = match_list_all
        d_mapping_resolved[id_protein_ref] = list(set(match_list_all))
        
        
        #== DEBUG:
        for i, (id_protein, m, mresolved) in enumerate(match_list_none):
            print 'match_list_none', i, id_protein, m, mresolved
        
        #== Write out match_list_none:
        self.output['match_list_none'] = match_list_none
            
            
        return d_mapping_resolved
    
    def get_d_mapping(self, fname_peptide_list_fasta, sim_cutoff_fraction, mappingAlgo):
        '''
        INPUT: (1) A list of peptides; (2) A list of polyproteins.
        OUTPUT: a dictionary containing for each protein, mapping results.
        '''
        d_mapping = mappingAlgo.search(fname_peptide_list_fasta, FNAME_PROTEINS, sim_cutoff_fraction = sim_cutoff_fraction) 
        #print 'debug: fname_peptide_list_fasta',  fname_peptide_list_fasta
        #print 'debug: FNAME_PROTEINS', FNAME_PROTEINS
        #print 'debug: d_mapping', d_mapping
        return d_mapping
    
    def get_d_mapping_sample(self):
        id_protein_list = ['0','1']
        d_mapping = {}
        match_list = [(85.0,  '1266', 1920, 1942),
                      (100.0, '1265', 1920, 1939),
                      (100.0, '1613', 396, 405),
                      (100.0, '82',   215, 228),
                      (100.0, '81',   215, 226)]
        d_mapping[id_protein_list[0]] = match_list
        
        match_list = [(90.0,  '1266', 1920, 1942),
                      (100.0, '1265', 1920, 1939),
                      (100.0, '1265', 100, 200),
                      (75,    '1613', 396, 405)]
        d_mapping[id_protein_list[1]] = match_list    
        return d_mapping

    def combine_mapping_list(self, d_mapping):
        '''
        d_id_protein[id_protein] = d_output_mapping
        d_output_mapping['mapping'][id_protein='']
        '''
        id_protein_list = d_mapping.keys(); id_protein_list.sort()
        match_list_all = []   # (id_protein, (sim, id_peptide, peptide, pos_i, pos_j))
        for id_protein in id_protein_list:
            match_list = d_mapping[id_protein]
            #print 'debug: match_list', match_list
            for m in match_list:
                row = (id_protein, m)
                match_list_all.append(row)
        match_list_all = sorted(match_list_all, key = lambda row: row[1][0], reverse=True)  # Sort from high to low similarity.
        return match_list_all

    def group_by_id_peptide(self, match_list_all):
        ''''''
        d_id_peptide = {}
        for match in match_list_all:
            #print 'match', match
            #(id_protein, (sim, id_peptide, peptide, pos_i, pos_j)) = match
            (id_protein, (sim, id_peptide, pos_i, pos_j)) = match
            d_id_peptide[id_peptide] = d_id_peptide.get(id_peptide, []) + [match]
            d_id_peptide[id_peptide] = sorted(d_id_peptide[id_peptide], key = lambda row: row[1][0], reverse=True)
        return d_id_peptide
    
    def group_by_id_protein(self, match_list_all):
        '''
        '''
        d_mapping = {}  # d[id_protein] = a list of 'row'; 'row' =  [sim_percent, id_peptide, pos_protein_i, pos_protein_j]
        for (id_protein, row) in match_list_all:
            d_mapping[id_protein] = d_mapping.get(id_protein, []) + [row]
        return d_mapping
    
    def filter_mapping_unique_keep_ties(self, d_id_peptide):
        '''
        (1) For each peptide, get a list of matches and sort from high similarity to low.
        (2) If highest match is unique, keep.
        (3) if highest match has a duplicate, then keep **only** duplicates.
        '''
        id_peptide_list = d_id_peptide.keys()
        for id_peptide in id_peptide_list:
            row_list = d_id_peptide[id_peptide]
            row_list_f = [row_list[0]]
            for i in range(1, len(row_list)):
                sim_before = row_list[i-1][1][0]   # (id_antigen, similarity, id_peptide, target_i, target_j) = row
                row = row_list[i]
                sim_current = row[1][0]
                if (sim_before > sim_current): break
                else: row_list_f.append(row)
            d_id_peptide[id_peptide] = row_list_f
        return d_id_peptide


def test_msa():
    '''Test reading of multiple sequence alignment using Biopython package'''
    msa = read_msa(FNAME_MSA)
    print msa
    for record in msa:
        print record.seq, '*', record.id

def test_pos_mapping_with_msa():
    '''Test a routine that returns positions with respect to msa.
    Appears to be working as expected.'''
    #==
    #pmap = PositionMappingWithMSA()
    #pmap.test_get_position_msa()
    
    #id_protein_ref = '1'
    #id_protein = '2'
    #pos_i = 3
    #pos_j = 5
    
    #id_protein_ref = '3'
    #id_protein = '2'
    #pos_i = 3
    #pos_j = 5
    
    #== This one maps to a gap in the reference protein. Should return 'None'
    id_protein_ref = '1'
    id_protein = '2'
    pos_i = 2
    pos_j = 5
      
    pmap = PositionMappingWithMSA()
    (pos_i_ref, pos_j_ref) = pmap.to_protein_reference(id_protein_ref, id_protein, pos_i, pos_j)
    print [id_protein,     pos_i,     pos_j]
    print [id_protein_ref, pos_i_ref, pos_j_ref]

def test_pos_mapping_with_msa_b():
    '''Check speed of two routines that return a dictionary for mapping between positions in sequences and MSA.'''
    import time
    
    seq = 'ABCDEF'
    seq_msa = 'ABCDEF'
    pmap = PositionMappingWithMSA()
    t1 = time.time()
    (da, db) = pmap.get_d_pos_seq_to_msa(seq, seq_msa)
    t2 = time.time()
    tdiff = t2 - t1
    print tdiff, (da,db)
    
    t1 = time.time()
    (da, db) = pmap.get_d_pos_seq_to_msa_fast(seq, seq_msa)
    t2 = time.time()
    tdiff = t2 - t1
    print tdiff, (da,db)

def test_msa_c():
    '''
    For a lsit of problematic peptides, map them to msa.
    '''
    fname_peptide_list = '/home/yohan/workspace/immunome_browser/experiments/20111129_hcv_mapping_stats/peptide_list_problematic.fasta'
    pmap = MappingWithMSA()
    d_mapping_resolved = pmap.search(fname_peptide_list, '', sim_cutoff_fraction=0.30)
    id_antigen_list = d_mapping_resolved.keys(); id_antigen_list.sort()
    for id_antigen in id_antigen_list:
        mlist = d_mapping_resolved[id_antigen]
        mlist = sorted(mlist, key=lambda row:row[0], reverse=True)
        for i, m in enumerate(mlist):
            print m     
   

if __name__ == '__main__':
    test_msa_c()
    #test_pos_mapping_with_msa_b()
    
    #test_pos_mapping_with_msa()
    #test_msa()
    
    #results  = test_run()
    #cPickle.dump(results, open('results.cpickle','w'))
    #(d_id_protein, match_list_all, d_id_peptide, match_list_all_f, d_output_mapping) = results
    