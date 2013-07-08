#! /usr/bin/python

import cPickle

'''There are a different classes of routines that I write very often:
    [] Reading and Writing routines.
    [] Checking routines.
    [] Analysis routines.
    []
    These boilerplate routines should be made into utilities.
    I tend to write too many variations of the same thing. Should minimize this.'''

#== Miscellaneous routines.
def remove_self(x, x_list):
    '''If x in in x_list, remove it. '''
    x_list_exclude = [xi for xi in x_list if xi != x]
    return x_list_exclude


def get_line_str(x, delimiter='\t'):
    x_str = ''
    for xi in x:
        x_str = x_str + str(xi)+delimiter
    return x_str

def get_str(v_list, delimiter=','):
    '''Returns a string version of a list surrounded by parentheses and separated by given delimiter.
       To be used in SQL query.'''
    v_list_str = ''
    for v in v_list:
        v_list_str = v_list_str + str(v) + delimiter   # Don't think there is a need to quote around epitope_id
    v_list_str = '(' + v_list_str[0:-1] + ')'
    return v_list_str

def write_cpickle(data, fname):
    f=open(fname,'w')
    cPickle.dump(data, f)
    f.close()

def read_cpickle(fname):
    f=open(fname,'r')
    data = cPickle.load(f)
    f.close()
    return data

def print_content_epitope_info(content):
    content_b = []
    for row in content:
        (epitope_id, epitope, source_accession, source_name) = row
        row_b = (int(source_accession), source_name, int(epitope_id), epitope)
        content_b.append(row_b)
    content_b.sort()
    for row in content_b:
        print row

def print_content(content):
    for line in content:
        print line

def write_file(content, fname):
    f = open(fname,'w')
    for row in content:
        f.write(row+'\n')
    f.close()

def get_content(fname):
    f=open(fname,'r')
    lines=f.readlines()
    f.close()
    return lines

def get_content_raw(fname):
    f=open(fname,'r')
    content=f.read()
    f.close()
    return content


