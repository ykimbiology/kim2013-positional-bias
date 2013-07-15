
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



#Example queries used in the project.
#All queries shown here were ran against database='iedb_public'. Last checked on July 1st 2013.
#For one source organism (i.e. vaccinia virus), the following sql queries demonstrate how to get different types of epitope data.


#source-organism='vaccinia virus'
SELECT eo.epitope_id, t.tcell_id, t.reference_id, t.as_char_value, t.as_num_responded, t.as_num_subjects
FROM tcell t, curated_epitope ce, epitope_object eo, organism o
WHERE t.curated_epitope_id = ce.curated_epitope_id
AND ce.e_object_id = eo.object_id
AND o.organism_id = eo.source_organism_org_id
AND (o.organism_id = '10245' OR o.path like '1:10239:35237:10240:10241:10242:10245:%')


#mhc-class='I'
SELECT eo.epitope_id, t.tcell_id, t.reference_id, t.as_char_value, t.as_num_responded, t.as_num_subjects
FROM tcell t, curated_epitope ce, epitope_object eo, organism o, mhc_allele_restriction m
WHERE t.curated_epitope_id = ce.curated_epitope_id
AND ce.e_object_id = eo.object_id
AND o.organism_id = eo.source_organism_org_id
AND (o.organism_id = '10245' OR o.path like '1:10239:35237:10240:10241:10242:10245:%')
AND m.displayed_restriction = t.mhc_allele_name
AND m.class = 'I'


# immunization-context='Organism'
SELECT eo.epitope_id, t.tcell_id, t.reference_id, t.as_char_value, t.as_num_responded, t.as_num_subjects
FROM tcell t, curated_epitope ce, epitope_object eo, organism o, mhc_allele_restriction m, object obj
WHERE t.curated_epitope_id = ce.curated_epitope_id
AND ce.e_object_id = eo.object_id
AND o.organism_id = eo.source_organism_org_id
AND (o.organism_id = '10245' OR o.path like '1:10239:35237:10240:10241:10242:10245:%')
AND m.displayed_restriction = t.mhc_allele_name
AND m.class = 'I'
AND t.iv1_imm_object_id=obj.object_id
AND obj.object_type ='Organism'

# immunization-context='peptide'
SELECT eo.epitope_id, t.tcell_id, t.reference_id, t.as_char_value, t.as_num_responded, t.as_num_subjects
FROM tcell t, curated_epitope ce, epitope_object eo, organism o, mhc_allele_restriction m, object obj
WHERE t.curated_epitope_id = ce.curated_epitope_id
AND ce.e_object_id = eo.object_id
AND o.organism_id = eo.source_organism_org_id
AND (o.organism_id = '10245' OR o.path like '1:10239:35237:10240:10241:10242:10245:%')
AND m.displayed_restriction = t.mhc_allele_name
AND m.class = 'I'
AND t.iv1_imm_object_id=obj.object_id
AND (obj.object_type='Fragment of a Natural Sequence Molecule' OR obj.object_type='Sequence Molecule No Natural Source')


# Given a list of epitope ids, how to retrieve source antigen information (e.g. source antigen ids, sequences):
# This sql returns source antigen information only for epitopes listed in (133, 570).
SELECT eo.epitope_id, obj.starting_position, obj.ending_position, s.source_id, s.accession, s.name, s.organism_name, s.sequence
FROM iedb_public.epitope_object eo, iedb_public.source s, iedb_public.object obj
WHERE eo.epitope_id in (133, 570)
AND eo.source_antigen_accession=s.accession
AND obj.object_id=eo.object_id
