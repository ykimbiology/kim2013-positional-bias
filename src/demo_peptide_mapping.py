#! /usr/bin/python -B

#Copyright (c) <2013>, <Yohan Kim>
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the <organization> nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
