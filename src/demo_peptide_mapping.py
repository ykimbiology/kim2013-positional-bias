#! /usr/bin/python

"""
Should run this script in the folder where this script is located.
"""

import sys
import os

from main.settings import *
from epitopemapping.mappeptides import map_peptides_selector



fname_peptides = '../examples/data/sequences/vacv/peptide_list.10.fasta'
fname_antigens = '../examples/data/sequences/vacv/sequences.10245.fasta'



#Using bruteforce approach to mapping peptides to antigens:
d_mapping, mappingAlgo = map_peptides_selector(fname_peptides, fname_antigens, method='bruteforce', sim_cutoff_fraction=1.0)

#Using blast:
#d_mapping, mappingAlgo = map_peptides_selector(fname_peptides, fname_antigens, method='blast_a', sim_cutoff_fraction=1.0)

id_antigen_list = d_mapping.keys()
print 'Number of antigens with mapped peptides', len(d_mapping)
print 'A list of mapped peptides for a given antigen', d_mapping[id_antigen_list[0]]

# Each entry contains (percent_seq_id, id_peptide, pos_begin, pos_end)
# ex: [(100.0, '232', 3, 11)]
