import time
import sys
import os
from subprocess import Popen, call


# For writing temporary files:
PATH_SCRATCH     = '/home/yohan/data_large/immunomebrowser/scratch'

# Where static html files will be written to. The django web will not use these.:
PATH_STATICFILES = os.path.join(PATH_SCRATCH, 'html')

# Where reference proteomes are located.
PATH_PROTEOME    = '/home/yohan/data_large/immunomebrowser/data/genomes/complete'



#== Specify location of the blast binaries.
PATH_BLAST   = '/home/yohan/resource_biology/blast/blast-2.2.23/bin'
CMD_BLAST_FORMATDB = os.path.join(PATH_BLAST, 'formatdb')
CMD_BLAST_BLASTALL = os.path.join(PATH_BLAST, 'blastall')

#== Specify location of the usearch binary:
# Example command: usearch6.0.152_i86linux32  -search_global peptide_list.10k.fasta -db sequences.10245.fasta  -id 0.5 -fulldp -alnout output.txt
PATH_USEARCH = '/home/yohan/resource_biology/usearch/usearch6.0'
CMD_USEARCH = os.path.join(PATH_USEARCH, 'usearch6.0.152_i86linux32')

#== smith-waterman SSEARCH
# Example command: ssearch36 -p peptide_list.fasta sequences.11103.fasta -m 8 -f 100 > output2.txt
PATH_SSEARCH = '/home/yohan/resource_biology/fasta/fasta-36.3.5a/bin'
CMD_SSEARCH = os.path.join(PATH_SSEARCH, 'ssearch36')

#== Specify db from which epitopes will be retrieved:
DB_EPITOPES = 'iedb_public'  # 'iedb_public', 'iedb_query' # 'iedb_query' used when curation-db specific feature is needed?


#== For plotting response frequencies == 
RES_NO_PLOT = 'NA'  # 'Residue for no plotting'. To be used by plot_factory in R. 'NA' will not be plotted.



#== Settings for peptide:antigen mapping:
# When running a commandline version, these params need to be set.
# If using the webapp, then the form allows the user to do this.
PEPTIDE_ANTIGEN_MAPPING_METHOD             = 'search_local_cmplt'    # [search_local_cmplt, search_local, blast_a, blast_b, bruteforce, hcv_msa]
PEPTIDE_ANTIGEN_SIMILARITY_CUTOFF_FRACTION = 0.5         # Default similarity cutoff of 80% sequence identity; For hcv, 50% was used.
ENFORCE_ANTIGEN_HOMOLOGY = False                          # When mapping epitopes to antigens, should source antigen be homologous to target antigens?
FILTER_MAPPING = False                                     # Should only best peptide:protein matches be kept? Or allow all peptide:protein mapping meeting the threshold? Ties are kept. [all, best_keep_ties, best_exclude_ties]





# To be consistent, use lowercase for the variables derived from the constants?
class SettingsMapping(object):
    '''Capture all mapping related settings information.'''
    def __init__(self, id_session = 'output'):
        self.peptide_antigen_mapping_method             = PEPTIDE_ANTIGEN_MAPPING_METHOD
        
        # ToDo: sim threshold are duplicated: Use only one.
        self.peptide_antigen_similarity_cutoff_fraction = PEPTIDE_ANTIGEN_SIMILARITY_CUTOFF_FRACTION  # [0,1]
        self.peptide_antigen_similarity_cutoff_percent  = 100.0*self.peptide_antigen_similarity_cutoff_fraction  # [0,100]
        self.enforce_antigen_homology                   = ENFORCE_ANTIGEN_HOMOLOGY
        self.filter_mapping                             = FILTER_MAPPING

        #== BLAST specific settings:
        self.cmd_formatdb                     = CMD_BLAST_FORMATDB
        self.cmd_blastall                     = CMD_BLAST_BLASTALL
        self.cmd_blast_peptides = ''  # blast command specific for peptides.
        self.cmd_blast_proteins = ''  # blast command specific for proteins.

        #== Location of complete proteomes of reference genomes:
        self.path_proteome                  = PATH_PROTEOME  # Used by util_antigens.py
        
        #== Where template html files are located.
        # To be used with django template; used in SQL query generation.
        self.path_templates                    = PATH_TEMPLATES
        
        #== Location of output files (html, csv, mapping data, etc.):
        self.path_staticfiles                  = PATH_STATICFILES # Those files for displaying on the web.

        '''#== The following directories will be created.
           Q: Should path_session be updated?'''
        self.id_session                       = id_session   # Without the directory prefix.
        self.path_session                     = os.path.join(self.path_staticfiles, id_session)                  # One Immunomebrowser session for one user.
        if DEBUG==False: self.path_session    = self.path_session + '_'+ str(time.time())  # To get unique directory in production.
        
        print 'debug path_session', self.path_session
        self.path_settings          = os.path.join(self.path_session,'settings')           # Store query-related information.
        self.path_peptides          = os.path.join(self.path_session,'peptides')
        self.path_antigens_homology = os.path.join(self.path_session,'antigens_homology')
        self.path_antigens_source   = os.path.join(self.path_session,'antigens_source')
        self.path_antigens_target   = os.path.join(self.path_session,'antigens_target')

        self.path_mapping           = os.path.join(self.path_session,'mapping')
        self.path_plotting          = os.path.join(self.path_session,'plotting')
        self.path_html              = os.path.join(self.path_session,'html')        

        self.clear_data()  # Removes the base output directory so that previous results do not accumulate.
        self.create_dirs()

    def get_dir_list(self):
        '''Make these directories if they are not present.'''
        dir_list = []
        dir_list.append(self.path_session)
        dir_list.append(self.path_settings)
        dir_list.append(self.path_peptides)
        dir_list.append(self.path_antigens_homology)
        dir_list.append(self.path_antigens_source)
        dir_list.append(self.path_antigens_target)
        dir_list.append(self.path_mapping)
        dir_list.append(self.path_plotting)
        dir_list.append(self.path_html)
        return dir_list

    def create_dirs(self):
        dir_list = self.get_dir_list()
        for dir_path in dir_list:
            if os.access(dir_path, os.F_OK) == False:  # I am not sure what 'mode' should be.
                os.makedirs(dir_path)
                Popen(['chmod', '777' ,'-R' ,dir_path])
                
    def clear_data(self, debug=False):
        ''' For a given list of directories, removes data within each.
            For now remove those in the output folder.
        '''
        if os.access(self.path_session, os.F_OK) == True:  # Not sure why this is returning True, eventhough the path does not exist. July 27 2010
            cmd = 'rm -fR ' + self.path_session
            if debug==True: print 'clear_data: ', self.path_session
            if debug==True: print 'clear_data: ', os.access(self.path_session, os.F_OK)
            if debug==True: print 'clear_data: ', cmd
            #p = Popen(cmd, shell=True)
            #p = Popen(cmd, shell=False) # This does not work.
            #p = call(cmd) # Does not work.
            p = call(['rm','-fR', self.path_session]) # This version works. Should use this because call waits for the process to end.
        else:
            if debug==True: print 'clear_data: Directory %s not found ' %(self.path_session)

    def __str__(self):
        entry_list = []
        entry_list.append('peptide_antigen_similarity_cutoff_fraction' +'\t' + str(self.peptide_antigen_similarity_cutoff_fraction))
        entry_list.append('peptide_antigen_mapping_method' +'\t' + str(self.peptide_antigen_mapping_method))
        entry_list.append('enforce_antigen_homology' +'\t' + str(self.enforce_antigen_homology))
        entry_list.append('filter_mapping' +'\t' + str(self.filter_mapping))

        line = '\n'.join(map(str, entry_list))
        return line



#== DEPRECATED ================================================================



if __name__ == '__main__':
    settings = SettingsMapping()
    settings.clear_data()