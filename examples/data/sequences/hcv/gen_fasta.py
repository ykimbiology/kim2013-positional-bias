#! /usr/bin/python

import csv
spamReader = csv.reader(open('table_export.csv', 'r'), delimiter=',', quotechar='\"')
for row in spamReader:
    id_peptide = row[0].strip()
    peptide = row[2].strip()
    if len(peptide)>0:
#        print [id_peptide, peptide]
        print '>id|'+id_peptide
        print peptide
   
