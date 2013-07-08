

'''Using cython to speed up for-loop.
'''
import numpy

def get_hamming_fraction_cython(peptide_a, peptide_b):
    '''Fraction of residues shared between the two peptides.
       One way to speed up is to use dot product of vectors? and sum 1s? Not sure how much faster?'''
    cdef int count_shared, peptide_len, i
    cdef float fraction_shared
    count_shared = 0
    peptide_len = len(peptide_a)
    for i in range(peptide_len):
        if (peptide_a[i] == peptide_b[i]): count_shared = count_shared + 1
    fraction_shared = float(count_shared)/peptide_len
    return fraction_shared


def get_hamming_fraction_numpy(peptide_a, peptide_b, peptide_len):
    '''Each peptide is a binary vector using numpy.array() '''
    count_match = numpy.dot(peptide_a, peptide_b)
    f_shared = float(count_match)/peptide_len
    return f_shared


def get_binary_vector(peptide_str):
    ''''''
    cdef int peptide_len, i, index
    aa_list = 'ACDEFGHIKLMNPQRSTVWY'
    peptide_len = len(peptide_str)
    vlen = 20*peptide_len
    v=numpy.zeros(vlen, int)
    for i,residue in enumerate(peptide_str):
        index = aa_list.index(residue)
        v[i*20+index] = 1
    return v

def get_match_list_cython(peptide_list, protein_seq, double sim_fraction_cutoff):
    '''HEre, id_peptide = index for each peptide. '''
    cdef int index_peptide, peptide_len, num_peptide, i, protein_seq_len
    cdef float f_score
    
    protein_seq_v = get_binary_vector(protein_seq)
    protein_seq_len = len(protein_seq)
    match_list = []  # row = (sim_percent, id_peptide, pos_protein_i, pos_protein_j)
    for index_peptide, peptide in enumerate(peptide_list):
        peptide_len = len(peptide)
        num_peptide = protein_seq_len - peptide_len + 1
        for i in range(num_peptide):
            peptide_test = protein_seq[i:i+peptide_len]
            f_score = get_hamming_fraction_cython(peptide, peptide_test)
            if f_score >= sim_fraction_cutoff:
                row = (100.0*f_score, index_peptide, i+1, i+peptide_len+1)
                match_list.append(row)
    return match_list



def get_match_list_numpy(peptide_list, protein_seq, double sim_fraction_cutoff):
    '''HEre, id_peptide = index for each peptide. '''
    cdef int index_peptide, peptide_len, num_peptide, i, protein_seq_len
    cdef float f_score
    
    protein_seq_v = get_binary_vector(protein_seq)
    protein_seq_len = len(protein_seq)
    match_list = []  # row = (sim_percent, id_peptide, pos_protein_i, pos_protein_j)
    for index_peptide, peptide in enumerate(peptide_list):
        peptide_v = get_binary_vector(peptide)
        peptide_len = len(peptide)
        num_peptide = protein_seq_len - peptide_len + 1
        for i in range(num_peptide):
            peptide_test = protein_seq[i:i+peptide_len]
            #peptide_test_v = get_binary_vector(peptide_test)
            peptide_test_v = protein_seq_v[i*20:i*20+peptide_len*20]
            f_score = get_hamming_fraction_numpy(peptide_v, peptide_test_v, peptide_len)
            if f_score >= sim_fraction_cutoff:
                row = (100.0*f_score, index_peptide, i+1, i+peptide_len+1)
                match_list.append(row)
    return match_list
