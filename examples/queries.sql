
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

