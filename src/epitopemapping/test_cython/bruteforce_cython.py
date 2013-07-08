
'''For testing purpose, brute force approach to peptide:protein mapping.'''
import random
import time

from Bio import SeqIO

def get_peptides(seq_antigen, peptide_len):
    peptide_list = []
    num_peptides = len(seq_antigen) - peptide_len + 1
    for i in range(num_peptides):
        peptide = seq_antigen[i:i + peptide_len]
        peptide_list.append(peptide)
    return peptide_list

def hamming_similarity(s1, s2):
    #assert len(s1) == len(s2)
    return sum(ch1 == ch2 for ch1, ch2 in zip(s1, s2))

def get_hamming_fraction(peptide_a, peptide_b):
    '''Fraction of residues shared between the two peptides.
       One way to speed up is to use dot product of vectors? and sum 1s? Not sure how much faster?'''
    count_shared = 0
    peptide_len = len(peptide_a)
    for i in range(peptide_len):
        res_a = peptide_a[i]
        res_b = peptide_b[i]
        if (res_a == res_b): count_shared = count_shared + 1
    fraction_shared = float(count_shared)/peptide_len
    return fraction_shared
    #== For some reason, the following is slower than above: 0.039 vs. 0.024
#    count_shared = hamming_similarity(peptide_a, peptide_b)
#    fraction_shared = float(count_shared)/len(peptide_a)
#    return fraction_shared


def get_peptide_list_random(n=1000, peptide_length=9):
    aa_list = 'ACDEFGHIKLMNPQRSTVWY'
    peptide_list = []
    #choose_residue_random = lambda
    for ni in range(n):
        peptide = [aa_list[random.randint(0,19)] for i in range(peptide_length)]
        peptide = ''.join(peptide)
        peptide_list.append(peptide)
    return peptide_list  

def get_protein_seq_random(protein_length=100):
    aa_list = 'ACDEFGHIKLMNPQRSTVWY'
    seq = ''.join([aa_list[random.randint(0,19)] for i in range(protein_length)])
    return seq
    

def get_match_list(peptide_list, protein_seq, sim_fraction_cutoff=1.0):
    '''HEre, id_peptide = index for each peptide. '''
    match_list = []  # row = (sim_percent, id_peptide, pos_protein_i, pos_protein_j)
    for index_peptide, peptide in enumerate(peptide_list):
        peptide_len = len(peptide)
        num_peptide = len(protein_seq) - peptide_len + 1
        for i in range(num_peptide):
            peptide_test = protein_seq[i:i+peptide_len]
            f_score = get_hamming_fraction(peptide, peptide_test)
            if f_score >= sim_fraction_cutoff:
                row = (100*f_score, index_peptide, i+1, i+peptide_len+1)
                match_list.append(row)
    return match_list
        
      

def get_hamming_fraction_fast(peptide_a, peptide_b, sim_fraction_cutoff):
    ''' While going through the residues, if current count shared is worse than the cutoff, break executed.
       Fraction of residues shared between the two peptides.
       One way to speed up is to use dot product of vectors? and sum 1s? Not sure how much faster?'''
    count_shared = 0
    count_missed = 0
    count_cutoff = sim_fraction_cutoff*len(peptide_a)
    for i in range(len(peptide_a)):
        res_a = peptide_a[i]
        res_b = peptide_b[i]
        if (res_a == res_b):
            count_shared = count_shared + 1
        else:
            count_missed = count_missed + 1
        if count_shared < count_cutoff: break
    fraction_shared = float(count_shared)/len(peptide_a)
    return fraction_shared

class BruteForce(object):
    '''INPUTS: (1) list of peptides (2) list of proteins.
       For each protein, generate overlapping windows of peptides and compare them against the input peptide.
       Repeat for all input peptides and return only those mappings that meet requirement.
       Returned residue positions start at 1, rather than 0.
       '''
    def __init__(self, debug=False):
        self.debug = debug
        self.output = None
    
    def read_fasta(self, fname):
        f = open(fname, 'r')
        record_list = list(SeqIO.parse(f, 'fasta'))
        f.close()
        return record_list
    
    def parse_ncbi_gi(self, record_protein):
        id_protein = record_protein.id.split('|')[1].strip()  # e.g. 'gi|1234'
        return id_protein
    
    def search(self, fname_peptides, fname_proteins, sim_cutoff_fraction=1.0):
        '''
        Loop over (peptides,proteins) and record each peptide:protein mapping that meets threshold.
            For each protein, generate appropriate set of peptides.
            
            row = (float(percent_identity), id_query, target_i_int, int(target_j))
            d_mapping[id_target] = d_mapping.get(id_target, []) + [row]
        '''
        record_list_peptide = self.read_fasta(fname_peptides)
        #== For debugging to speed up comparison:
        #record_list_peptide = record_list_peptide[0:100]
        
        record_list_protein = self.read_fasta(fname_proteins)
        d_mapping = self.search_records(record_list_peptide, record_list_protein, sim_cutoff_fraction=sim_cutoff_fraction)
        return d_mapping

    def search_records(self, record_list_peptide, record_list_protein, sim_cutoff_fraction=1.0):
        '''
        This version uses BioPythons SeqRecord object:
        Loop over (peptides,proteins) and record each peptide:protein mapping that meets threshold.
            For each protein, generate appropriate set of peptides.
            
            row = (float(percent_identity), id_query, target_i_int, int(target_j))
            d_mapping[id_target] = d_mapping.get(id_target, []) + [row]
        NOTE: Returned residues use 1-index. First residue gets a value of 1.
        '''
        d_mapping = {}
        
        for (i,record_protein) in enumerate(record_list_protein):
            seq_protein = record_protein.seq
            id_protein = self.parse_ncbi_gi(record_protein)
            
            #== Mapping a peptide to a sequence. This section should be optimized for speed:
            if self.debug==True: record_list_peptide = record_list_peptide[0:100]
            
            #== Looping over a list of peptides from iedb:
            # INPUT: (1) list of peptides, (2) protein sequence.
            # OUTPUT: A list of matches. row = (sim_percent, id_peptide, pos_protein_i, pos_protein_j)
            # get_match_list(peptide_list, protein_seq)
            for (j, record_peptide) in enumerate(record_list_peptide):  
                seq_peptide = record_peptide.seq
                id_peptide = self.parse_ncbi_gi(record_peptide)
                #print 'seq_peptide', seq_peptide
                if self.debug==True: 
                    if j % 1000 == 0: print ['search_records', i, j, len(record_list_peptide), id_peptide]
                temp_peptide_list = get_peptides(seq_protein, len(seq_peptide)) 
                sim_list_fraction = [get_hamming_fraction(temp_peptide, seq_peptide) for temp_peptide in temp_peptide_list]
                for pos in range(len(temp_peptide_list)):
                    temp_peptide = temp_peptide_list[pos]  # This peptide is from seq_protein, NOT the query peptide.
                    sim_fraction = sim_list_fraction[pos]
                    if (sim_fraction >= sim_cutoff_fraction):
                        sim_percent = 100.0*sim_fraction
                        pos_protein_i = pos + 1                  # Residue number starts at 1, rather than 0.
                        pos_protein_j = pos_protein_i + len(seq_peptide) - 1   # index_a = 1, index_b = 10, len(peptide)=10
                        row = (sim_percent, id_peptide, pos_protein_i, pos_protein_j)
                        d_mapping[id_protein] = d_mapping.get(id_protein, []) + [row]
            if self.debug==True: print '\t'.join(map(str,[i, id_protein, len(d_mapping.get(id_protein,[]) )]))
        return d_mapping

def test_speed_hamming():
    '''
    Find out which version of peptide similarity calculation is faster.
    '''
    n = 10000
    peptide_list_a = get_peptide_list_random(n=n)
    peptide_list_b = get_peptide_list_random(n=n)
    
    t1 = time.time()
    sim_list_a = []
    for (pi,pj) in zip(peptide_list_a, peptide_list_b):
        similarity = get_hamming_fraction(pi, pj)
        sim_list_a.append(similarity)
    t2 = time.time()
    tdiff_a = t2-t1
    
    
    t1 = time.time()
    sim_list_b = []
    sim_fraction_cutoff = 0.30
    for (pi,pj) in zip(peptide_list_a, peptide_list_b):
        similarity = get_hamming_fraction_fast(pi,pj,sim_fraction_cutoff)
        sim_list_b.append(similarity)
    t2 = time.time()
    tdiff_b = t2-t1

    print ['before', tdiff_a]
    print ['after', tdiff_b]
    
    
def test_mapping():
    '''
    Demonstrate mapping of peptides to proteins using BruteForce class.
    '''
    fname_peptides = ''   # A fasta file containing peptide sequences.
    fname_proteins = ''   # A fasta file containing protein sequences.
    bf = BruteForce()
    d_mapping = bf.search(fname_peptides, fname_proteins, sim_cutoff_fraction=1.0)
    id_protein_list = d_mapping.keys()
    print 'Number of proteins ', len(id_protein_list)

def test_hamming_speed_01():
    import time

    #Importing cython module:
    import pyximport; pyximport.install()
    import hamming_cython

    num_peptide = 1000
    protein_length = 1000
    peptide_list = get_peptide_list_random(n=num_peptide)
    protein_seq = get_protein_seq_random(protein_length=protein_length)

    xi = time.time()
    match_list = get_match_list(peptide_list, protein_seq, sim_fraction_cutoff=0.5)
    xj = time.time()
    dx = xj-xi
    print ['%.2fs' % dx, 'num_peptide=',num_peptide, 'protein_length=',protein_length, 'len(match_list)', len(match_list)]
    #print match_list[0:5]

    xi = time.time()
    match_list = hamming_cython.get_match_list_cython(peptide_list, protein_seq, 0.5)
    xj = time.time()
    dx = xj-xi
    print ['%.2fs' % dx, 'num_peptide=',num_peptide, 'protein_length=',protein_length, 'len(match_list)', len(match_list)]
    #print match_list[0:5]

    xi = time.time()
    match_list = hamming_cython.get_match_list_numpy(peptide_list, protein_seq, 0.5)
    xj = time.time()
    dx = xj-xi
    print ['%.2fs' % dx, 'num_peptide=',num_peptide, 'protein_length=',protein_length, 'len(match_list)', len(match_list)]
    #print match_list[0:5]

if __name__ == '__main__':
    #test_speed_hamming()
    test_hamming_speed_01()
