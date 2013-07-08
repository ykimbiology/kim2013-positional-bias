#! /usr/bin/python


fname='table_export.csv'



import csv
reader = csv.reader(open(fname,'r'), delimiter=',', quotechar="\"")
for row in reader:
    id_peptide = row[0].strip()
    peptide = row[1].strip()
    print '>id|'+id_peptide
    print peptide
