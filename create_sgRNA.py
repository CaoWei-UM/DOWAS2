#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 06:57:11 2018

@author: chenrui
"""

import itertools
import os
import sys


sgRNA = sys.argv[1]
output_file_location = sys.argv[2]
sequence = sgRNA
output_path = output_file_location

def write_library(output_path, library,mismatch,name):
    with open(output_path,'a') as f:
        for i,sequence in enumerate(library):
            f.write('>{0}_{1}{2}_{3}\n'.format(mismatch,name,i,i))
            f.write(sequence + '\n')

def clear_file(output_path):
    if os.path.exists(output_path):
        os.remove(output_path)

def neighborhood(word,mismatch):
    seed = []
    word_list = list(itertools.product('ATCG',repeat=mismatch))
    for site in itertools.combinations(range(len(word)),mismatch):
        for product in word_list:
            word_copy = list(word)
            flag = True
            for i in range(mismatch):
                if word_copy[site[i]] == product[i]:
                    flag = False
                    break
                else:
                    word_copy[site[i]] = product[i]
            if flag:
                seed.append(''.join(word_copy))
    return seed

def create_library(seed,mismatch,pam):
    library = []
    for i in 'ATCG':
        text = ['{0}{1}{2}'.format(k,i,pam) for k in seed]
        library.extend(text)
    return library

clear_file(output_path)
for i in range(1,7):
    seed = neighborhood(sequence,i)
    library_NGG = create_library(seed,i,'GG')
    write_library(output_path,library_NGG,i,'NGG')    
    if i < 5:
        library_NAG = create_library(seed,i,'AG')
        write_library(output_path,library_NAG,i,'NAG')
