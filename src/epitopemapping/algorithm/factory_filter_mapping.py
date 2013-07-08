
import sys
import time

'''
GOAL: Mapping records can be filtered using various rules.
'''


def group_by_id_peptide(d_mapping_antigen):
    d_mapping_peptide = {}
    id_target_list = d_mapping_antigen.keys()
    for id_target in id_target_list:
        mapped_peptide_list = d_mapping_antigen[id_target]
        for row in mapped_peptide_list:
            (similarity, id_peptide, target_i, target_j) = row
            row_b = (id_target, similarity, id_peptide, target_i, target_j)
            d_mapping_peptide[id_peptide] = d_mapping_peptide.get(id_peptide, []) + [row_b]
    id_peptide_list = d_mapping_peptide.keys()
    for id_peptide in id_peptide_list:
        d_mapping_peptide[id_peptide] = sorted(d_mapping_peptide[id_peptide], key = lambda row: float(row[1]), reverse=True)
    return d_mapping_peptide

def group_by_id_antigen(d_mapping_peptide):
    d_mapping_antigen = {}
    id_peptide_list = d_mapping_peptide.keys()
    for id_peptide in id_peptide_list:
        row_list = d_mapping_peptide[id_peptide]
        for row in row_list:
            (id_antigen, similarity, id_peptide, target_i, target_j) = row
            row_b = (similarity, id_peptide, target_i, target_j)
            d_mapping_antigen[id_antigen] = d_mapping_antigen.get(id_antigen, []) + [row_b]
    return d_mapping_antigen


# TODO
## Write a python decoraor to measure time a routine takes?
#ti = time.time()
#tj = time.time(); td = tj-ti; 
#report_time(ti,tj, '== Map peptides: filter_mapping_by_homology_protein ==')
#    
#
#def measure_time(fn):
#    """
#    Python decorator.
#    """
#    def wrapped():
#        return fn
#    return wrapped


def filter_mapping_by_homology_protein(d_mapping,
                                       d_peptide_info,
                                       fname_antigens_db_source, 
                                       fname_antigens_db_target):
    """
    Filter mapping results based on homology between source and target antigens.
    """
    from blast import HomologousObject
    hobj = HomologousObject(fname_antigens_db_source, fname_antigens_db_target)
    d_mapping_f = filter_mapping_with_homology(d_mapping, d_peptide_info, hobj)
    return (d_mapping_f, hobj)

def filter_mapping_with_homology(d_mapping, d_peptide_info, hobj):
    """ 
    GOAL: Exclude those peptide:protein mapping that involve nonhomologous source:target protein relationship based on the provided hobj.
        d_peptide_info: For each peptide, there is a list of id_antigen of source antigens. This list will be compared against the blast hits.
        
        Not sure whether to make a deep copy of the dictionary: d_mapping_results
        == Pseudo Programming Process ==
        Loop over each id_antigen_target:
            For each id_target, get a list of all mapped peptides.
            For each peptide, return a list of its id_antigens_source.
            Exclude peptide mapping if antigen_source and antigen_target are not homologous.
    """
    id_target_list = d_mapping.keys()
    
    for id_target in id_target_list:
        mapped_peptide_list = d_mapping[id_target]
        mapped_peptide_list_f = []  # Filtered based on antigen homology.
        for row in mapped_peptide_list:
            #== get a list of id_antigen_source for each peptide.
            #== Then determine whether any of it is homologous to id_antigen_target.
            (similarity, id_peptide, target_i, target_j) = row
            d_row_list = d_peptide_info[id_peptide]
            #id_source_accession_list = [d_row['source_accession'] for d_row in d_row_list]
#            print 'DEBUG: filter_mapping_with_homology len(d_row_list)', len(d_row_list), d_row_list
#            if len(d_row_list)>0:
#                print 'DEBUG: filter_mapping_with_homology d_row_list[0]', d_row_list[0]
            #id_source_list = [d_row['source_id'] for d_row in d_row_list]
            id_source_list = [d_row['source_accession'] for d_row in d_row_list] # Antigen-homology feature will use source_accession to identify source antigens.
            is_homologous = hobj.check_homologous_list(id_source_list, id_target)
            
            #print 'DEBUG', [is_homologous, id_target, id_source_list]
            if is_homologous == True: mapped_peptide_list_f.append(row)
        d_mapping[id_target] = mapped_peptide_list_f
    return d_mapping


def filter_mapping_unique_keep_ties(d_mapping_antigen):
    """
    (1) For each peptide, get a list of matches and sort from high similarity to low.
    (2) If highest match is unique, keep.
    (3) if highest match has a duplicate, then keep **only** duplicates.
    """
    d_mapping_peptide = group_by_id_peptide(d_mapping_antigen)
    id_peptide_list = d_mapping_peptide.keys()
    for id_peptide in id_peptide_list:
        row_list = d_mapping_peptide[id_peptide]
        row_list_f = [row_list[0]]
        for i in range(1, len(row_list)):
            sim_before = row_list[i-1][1]   # (id_antigen, similarity, id_peptide, target_i, target_j) = row
            row = row_list[i]
            sim_current = row[1]
            if (sim_before > sim_current): break
            else: row_list_f.append(row)
        d_mapping_peptide[id_peptide] = row_list_f
    d_mapping_anigen_b = group_by_id_antigen(d_mapping_peptide)
    return d_mapping_anigen_b

def filter_mapping_unique_exclude_ties(d_mapping_antigen):
    """
    Tested on the example d_mapping_antigen data.
    """
    d_mapping_peptide = group_by_id_peptide(d_mapping_antigen)
    id_peptide_list = d_mapping_peptide.keys()
    for id_peptide in id_peptide_list:
        row_list = d_mapping_peptide[id_peptide]
        row_list_f = []
        if len(row_list)==1:
            row_list_f = [row_list[0]]
        else:
            sim_before = row_list[0][1]   # (id_antigen, similarity, id_peptide, target_i, target_j) = row
            sim_current = row_list[1][1]
            if (sim_before > sim_current): # If first simimilarity is higher than the second, include the first.
                row_list_f = [row_list[0]]
            else: row_list_f = []
        print 'debug', len(row_list_f)
        if len(row_list_f) == 0:
            d_mapping_peptide.__delitem__(id_peptide)
        else:
            d_mapping_peptide[id_peptide] = row_list_f
    d_mapping_anigen_b = group_by_id_antigen(d_mapping_peptide)
    return d_mapping_anigen_b