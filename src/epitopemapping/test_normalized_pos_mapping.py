#! /usr/bin/python


'''
Q: After peptides are sampled from a set of proteins and mapped to proteins using blast,
    what is the distribution of their normalized position?
    If this is not uniform, then there is a serious flaw in the mapping approach:
    
    [] Given the list of peptide list retrieved from iedb, sample a set of peptides uniform from
       the set of proteins so that distribution of lengths is the same.
    
    list_peptide_sampled = []
    for peptide in list_peptide:
        peptide_length = len(peptide)
        peptide_sampled = sample_peptide_from_proteins(peptide_length=9, seq_protein):
        list_peptide_sampled.append(peptide_sampled)
    
    calculate_normalized_position(d_mapping):
        #GOAL: For each peptide mapping, return normalized position score.
        
    OUTPUT: (1) A list of sampled peptides.
            (2) The set of proteins used.
            (3) Mapping results
            (4) Using (3), calculate a distribution of normalized positions.
'''
import random
import cPickle
import os

from Bio import SeqIO

FNAME_LIST_PROTEIN = '/home/yohan/workspace/immunome_browser/src/immunomebrowser.epitopemapping/immunomebrowser/epitopemapping/test_data/antigens_target/sequences.1773.fasta'
FNAME_LIST_PEPTIDE = '/home/yohan/workspace/immunome_browser/src/immunomebrowser.epitopemapping/immunomebrowser/epitopemapping/test_data/peptides_info/peptide_list.fasta'

class ProteinSet(object):
    def __init__(self):
        self.seq_list = [seq_record for seq_record in SeqIO.parse(FNAME_LIST_PROTEIN, "fasta")]
        #print self.seq_list[0]
        self.d_seq = {}
        for seq in self.seq_list:
            id_seq = self.parse_ncbi_gi(seq.id)
            self.d_seq[id_seq] = str(seq.seq)
            

    def parse_ncbi_gi(self, seq_id):
        id_seq = seq_id.split('|')[1].strip()
        return id_seq
        
    #public
    def get_random(self):
        ''''''
        seq_record = random.choice(self.seq_list)
        seq_protein = seq_record.seq
        #print seq_protein
        return str(seq_protein)
    
    def get_length(self, id_seq):
        return len(self.d_seq[id_seq])

class PeptideSet(object):
    def __init__(self):
        self.seq_list = [seq_record for seq_record in SeqIO.parse(FNAME_LIST_PEPTIDE, "fasta")]
        #print self.seq_list[0]
        self.d_seq = {}
        for seq in self.seq_list:
            id_seq = self.parse_ncbi_gi(seq.id)
            self.d_seq[id_seq] = str(seq.seq)
            
    def get_list_length(self):
        ''''''
        peptide_length_list = [len(seq) for seq in self.seq_list]
        return peptide_length_list

    def parse_ncbi_gi(self, seq_id):
        id_seq = seq_id.split('|')[1].strip()
        return id_seq
    
    def get_length(self, id_seq):
        return len(self.d_seq[id_seq])

def sample_peptide_from_proteins(seq_set, peptide_length=9):
    '''
    '''
    seq_protein = seq_set.get_random()
    seq_protein_len = len(seq_protein)
    index_range = range(seq_protein_len - peptide_length + 1)  # ex: for 9mer, range(9 - 9 + 1) = [0]
    index_peptide_random = random.choice(index_range)
    peptide_random = seq_protein[index_peptide_random: index_peptide_random + peptide_length]
    return peptide_random
    
def write_fasta(seq_list, fname_output):
    f=open(fname_output, 'w')
    for i,seq in enumerate(seq_list):
        f.write('>id|'+str(i+1)+'\n')
        f.write(seq+'\n')
    f.close()

def generate_peptide_list_sampled():
    peptide_set = PeptideSet()
    protein_set = ProteinSet()
    plen_list = peptide_set.get_list_length()
    print plen_list[0:10]
    peptide_list_sampled = []
    for plen in plen_list:
        peptide = sample_peptide_from_proteins(protein_set, peptide_length=plen)
        #print peptide
        peptide_list_sampled.append(peptide)
    print len(peptide_list_sampled), peptide_list_sampled[0:5]  
    write_fasta(peptide_list_sampled, 'peptide_list_sampled.fasta')
    
    
def get_fraction_position_bp_b(pos_start, pos_end, seq_length):
    '''Second version suggested by bpeters.
       pos_start uses zero-index system.'''
    peptide_length = pos_end - pos_start + 1
    pos_fraction = float(pos_start)/(seq_length - peptide_length + 1)
    assert (pos_fraction <= 1.0) and (pos_fraction >= 0.0)
        #print 'error pos fraction', pos_start, pos_end, seq_length, pos_fraction
    return pos_fraction

def get_normpos_list(fname_mapping):
    protein_set = ProteinSet()
    
    #fname_mapping = 'd_mapping.blast_a.sampled.cpickle'
    #fname_mapping = 'd_mapping.bruteforce.sampled.cpickle'
    d_mapping = cPickle.load(open(fname_mapping,'r'))
    id_protein_list = d_mapping.keys()
    content = []
    header = ['normpos', 'peptide_len', 'protein_len']
    content.append(header)
    for id_protein in id_protein_list:
        mapping_list = d_mapping[id_protein]
        protein_len = protein_set.get_length(id_protein)
        for (sim, id_peptide, pos_start, pos_end) in mapping_list:
            peptide_len = pos_end - pos_start + 1
            pos_normalized = get_fraction_position_bp_b(pos_start, pos_end, protein_len)
            print '\t'.join(map(str,[pos_normalized, peptide_len, protein_len]))
            row = map(str,[pos_normalized, peptide_len, protein_len])
            content.append(row)
    return content

def generate_normpos_set(method = 'blast_a', num_set=1):
    '''
    '''
    from test_map_peptides import generate_d_mapping
    from util_blastpeptides import write_file
    
    # sample peptides.
    # map peptides.
    # write normpos data
    # repeat N times.
    for i in range(num_set):
        fname_peptide_list = 'peptide_list_sampled.fasta'
        fname_mapping = '.'.join(['d_mapping',method,'cpickle'])
        generate_peptide_list_sampled()
        generate_d_mapping(method = method, fname_peptide_list=fname_peptide_list)
        content = get_normpos_list(fname_mapping)
        fname_normpos = '.'.join(['normpos', method, str(i), 'txt'])
        write_file(['\t'.join(row) for row in content], fname_normpos)
        os.remove(fname_peptide_list)
        os.remove(fname_mapping)
    
generate_normpos_set(method = 'bruteforce', num_set=10)
generate_normpos_set(method = 'blast_a', num_set=10)