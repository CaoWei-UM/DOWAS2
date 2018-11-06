# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import vcf
import pysam
import os
import sys

def write_haplotype(output_path, haplotype):
    with open(output_path,'a') as f:
        for i, sequence in haplotype.items():
            f.write('>' + str(i) + '\n')
            f.write(sequence + '\n')

def write_haplotype_name(name_output_path, haplotype_name):
    with open(name_output_path,'a') as f:
        for i, name in haplotype_name.items():
            f.write(str(i) + '>')
            f.write(name + '\n')

def clear_file(output_path):
    if os.path.exists(output_path):
        os.remove(output_path)

def construct_haplotype(reference, variants, sample, chromosome, start, end=None):
    '''
    construct haplotype at the special locate from the reference and sample variants.
    ''' 
    sequence = list(str(reference.fetch(chromosome, start, end)).upper())
    left_sequence = sequence.copy()
    right_sequence = sequence.copy()
    left_name = '>{0}:{1}:0_0'.format(chromosome,start)
    right_name = '>{0}:{1}:0_0'.format(chromosome,start)
    left_indel_length = 0
    right_indel_length = 0
    left_move = 0
    right_move = 0
    left_snp_name = ''
    right_snp_name = ''
    for record in variants.fetch(chromosome, start, end):
        call = record.genotype(sample)
        alleles = [str(i) for i in record.alleles]
        alleles_length = [len(i) for i in alleles]
        move = record.start - start

        if call.phased:
            (left, right) = [int(i) for i in call['GT'].split("|")]
#        elif 'PGT' in str(call):
#            (left, right) = [int(i) for i in call['PGT'].split("|")]
        else:
            (left, right) = [int(i) for i in call['GT'].split("/")]
        if alleles[0] == "".join(left_sequence[move:move+alleles_length[0]]):
            if left != 0:
                if alleles_length[0] > alleles_length[left]:
                    for i in range(alleles_length[0]):                    
                        left_sequence[move + i] = ''
                    for i in range(alleles_length[left]):                    
                        left_sequence[move + i] = alleles[left][i]
                elif alleles_length[0] < alleles_length[left]:
                    left_sequence[move] = alleles[left][0:alleles_length[left]-alleles_length[0]+1]
                else:
                    left_sequence[move] = alleles[left][0:1]
                left_indel_length = alleles_length[0] - alleles_length[left]
                if left_indel_length != 0:
                    left_move = left_indel_length                
                    left_name = '{0}+{1}_{2}'.format(left_name, move , left_move)
                else:
                    if left_snp_name:
                        left_snp_name = '{0}+{1}_{2}'.format(left_snp_name, move, alleles[left][0:1])
                    else:
                        left_snp_name = '{0}_{1}'.format(move, alleles[left][0:1])   
        if alleles[0] == "".join(right_sequence[move:move+alleles_length[0]]):
            if right != 0:
                if alleles_length[0] > alleles_length[right]:
                    for i in range(alleles_length[0]):                    
                        right_sequence[move + i] = ''
                    for i in range(alleles_length[right]):                    
                        right_sequence[move + i] = alleles[right][i]
                elif alleles_length[0] < alleles_length[right]:
                    right_sequence[move] = alleles[right][0:alleles_length[right]-alleles_length[0]+1]
                else:
                    right_sequence[move] = alleles[right][0:1]
                right_indel_length = alleles_length[0] - alleles_length[right]
                if right_indel_length != 0:
                    right_move = right_indel_length
                    right_name = '{0}+{1}_{2}'.format(right_name, move , right_move)                       
                else:
                    if right_snp_name:
                        right_snp_name = '{0}+{1}_{2}'.format(right_snp_name, move, alleles[right][0:1])
                    else:
                        right_snp_name = '{0}_{1}'.format(move, alleles[right][0:1])
    if left_snp_name:
        left_name = left_name + ':' + left_snp_name
    if right_snp_name:
        right_name = right_name + ':' + right_snp_name
    left_sequence = ''.join(left_sequence)
    right_sequence = ''.join(right_sequence)  
    return left_name, right_name, left_sequence, right_sequence

def construct_haplotype_state5(reference, variants, sample, chromosome, start, end=None):
    '''
    construct haplotype list at the special locate from the reference and sample variants.
    ''' 
    sequence = list(str(reference.fetch(chromosome, start, end)).upper())
    left_sequence = [sequence.copy()]
    right_sequence = [sequence.copy()]
    left_name = ['>{0}:{1}:0_0'.format(chromosome,start)]
    right_name = ['>{0}:{1}:0_0'.format(chromosome,start)]
    left_indel_length = 0
    right_indel_length = 0
    left_move = [0]
    right_move = [0]
    left_snp_name = ['']
    right_snp_name = ['']
    for record in variants.fetch(chromosome, start, end):
        call = record.genotype(sample)
        left_sequence_num = len(left_sequence)
        right_sequence_num = len(right_sequence)
        alleles = [str(i) for i in record.alleles]
        alleles_length = [len(i) for i in alleles]
        move = record.start - start
        
        if call.phased:
            (left, right) = [int(i) for i in call['GT'].split("|")]
#            print(move,' ',left,'|',right)
            combination = 1
        else:
            (left, right) = [int(i) for i in call['GT'].split("/")]   
            if left != right:
                combination = 2  
            else:
                combination = 1                
        if combination == 1:
            for j in range(left_sequence_num):
                if alleles[0] == "".join(left_sequence[j][move:move+alleles_length[0]]):
                    if left != 0:
                        if alleles_length[0] > alleles_length[left]:
                            for i in range(alleles_length[0]):                    
                                left_sequence[j][move + i] = ''
                            for i in range(alleles_length[left]):                    
                                left_sequence[j][move + i] = alleles[left][i]
                        elif alleles_length[0] < alleles_length[left]:
                            left_sequence[j][move] = alleles[left][0:alleles_length[left]-alleles_length[0]+1]
                        else:
                            left_sequence[j][move] = alleles[left][0:1]
                        left_indel_length = alleles_length[0] - alleles_length[left]
                        if left_indel_length != 0:
                            left_move[j] = left_indel_length                
                            left_name[j] = '{0}+{1}_{2}'.format(left_name[j], move , left_move[j])
                        else:
                            if left_snp_name[j]:
                                left_snp_name[j] = '{0}+{1}_{2}'.format(left_snp_name[j], move, alleles[left][0:1])
                            else:
                                left_snp_name[j] = '{0}_{1}'.format(move, alleles[left][0:1])
#                       if move == 56959:
#                           print(left_name)
            for j in range(right_sequence_num):
                if alleles[0] == "".join(right_sequence[j][move:move+alleles_length[0]]):    
                    if right != 0:
                        if alleles_length[0] > alleles_length[right]:
                            for i in range(alleles_length[0]):                    
                                right_sequence[j][move + i] = ''
                            for i in range(alleles_length[right]):                    
                                right_sequence[j][move + i] = alleles[right][i]
                        elif alleles_length[0] < alleles_length[right]:
                            right_sequence[j][move] = alleles[right][0:alleles_length[right]-alleles_length[0]+1]
                        else:
                            right_sequence[j][move] = alleles[right][0:1]
                        right_indel_length = alleles_length[0] - alleles_length[right]
                        if right_indel_length != 0:
                            right_move[j] = right_indel_length
                            right_name[j] = '{0}+{1}_{2}'.format(right_name[j], move , right_move[j])                       
                        else:
                            if right_snp_name[j]:
                                right_snp_name[j] = '{0}+{1}_{2}'.format(right_snp_name[j], move, alleles[right][0:1])
                            else:
                                right_snp_name[j] = '{0}_{1}'.format(move, alleles[right][0:1])
#                       if move == 12233:
#                           print(left_name)
#                           print(right_name)
        elif combination == 2:
#            print(left_name)
#            print(right_name)
            for j in range(left_sequence_num):                  
                left_sequence.append(left_sequence[j].copy())
                left_move.append(left_move[j])
                left_name.append(left_name[j])
                left_snp_name.append(left_snp_name[j])
                if alleles[0] == "".join(left_sequence[j][move:move+alleles_length[0]]):
                    if left != 0:
                        if alleles_length[0] > alleles_length[left]:
                            for i in range(alleles_length[0]):                    
                                left_sequence[j][move + i] = ''
                            for i in range(alleles_length[left]):                    
                                left_sequence[j][move + i] = alleles[left][i]
                        elif alleles_length[0] < alleles_length[left]:
                            left_sequence[j][move] = alleles[left][0:alleles_length[left]-alleles_length[0]+1]
                        else:
                            left_sequence[j][move] = alleles[left][0:1]
                        left_indel_length = alleles_length[0] - alleles_length[left]
                        if left_indel_length != 0:
                            left_move[j] = left_indel_length                
                            left_name[j] = '{0}+{1}_{2}'.format(left_name[j], move , left_move[j])
                        else:
                            if left_snp_name[j]:
                                left_snp_name[j] = '{0}+{1}_{2}'.format(left_snp_name[j], move, alleles[left][0:1])
                            else:
                                left_snp_name[j] = '{0}_{1}'.format(move, alleles[left][0:1])
                    if right != 0:
                        k = j + left_sequence_num
                        if alleles_length[0] > alleles_length[right]:
                            for i in range(alleles_length[0]):                    
                                left_sequence[k][move + i] = ''
                            for i in range(alleles_length[right]):                    
                                left_sequence[k][move + i] = alleles[right][i]
                        elif alleles_length[0] < alleles_length[right]:
                            left_sequence[k][move] = alleles[right][0:alleles_length[right]-alleles_length[0]+1]
                        else:    
                            left_sequence[k][move] = alleles[right][0:1]
                        right_indel_length = alleles_length[0] - alleles_length[right]
                        if right_indel_length != 0:
                            left_move[k] = right_indel_length
                            left_name[k] = '{0}+{1}_{2}'.format(left_name[k], move  , left_move[k])                       
                        else:
                            if left_snp_name[k]:
                                left_snp_name[k] = '{0}+{1}_{2}'.format(left_snp_name[k], move, alleles[right][0:1])
                            else:
                                left_snp_name[k] = '{0}_{1}'.format(move, alleles[right][0:1]) 
            for j in range(right_sequence_num):                  
                right_sequence.append(right_sequence[j].copy())
                right_move.append(right_move[j])
                right_name.append(right_name[j])
                right_snp_name.append(right_snp_name[j])
                if alleles[0] == "".join(right_sequence[j][move:move+alleles_length[0]]):
                    if right != 0:
                        if alleles_length[0] > alleles_length[right]:
                            for i in range(alleles_length[0]):                    
                                right_sequence[j][move + i] = ''
                            for i in range(alleles_length[right]):                    
                                right_sequence[j][move + i] = alleles[right][i]
                        elif alleles_length[0] < alleles_length[right]:
                            right_sequence[j][move] = alleles[right][0:alleles_length[right]-alleles_length[0]+1]
                        else:
                            right_sequence[j][move] = alleles[right][0:1]
                        right_indel_length = alleles_length[0] - alleles_length[right]
                        if right_indel_length != 0:
                            right_move[j] = right_indel_length
                            right_name[j] = '{0}+{1}_{2}'.format(right_name[j], move  , right_move[j])                       
                        else:
                            if right_snp_name[j]:
                                right_snp_name[j] = '{0}+{1}_{2}'.format(right_snp_name[j], move, alleles[right][0:1])
                            else:
                                right_snp_name[j] = '{0}_{1}'.format(move, alleles[right][0:1]) 
                    if left != 0:
                        k = j + right_sequence_num
                        if alleles_length[0] > alleles_length[left]:
                            for i in range(alleles_length[0]):                    
                                right_sequence[k][move + i] = ''
                            for i in range(alleles_length[left]):                    
                                right_sequence[k][move + i] = alleles[left][i]
                        elif alleles_length[0] < alleles_length[left]:
                            right_sequence[k][move] = alleles[left][0:alleles_length[left]-alleles_length[0]+1]
                        else:
                            right_sequence[k][move] = alleles[left][0:1]
                        left_indel_length = alleles_length[0] - alleles_length[left]
                        if left_indel_length != 0:
                            right_move[k] = left_indel_length                
                            right_name[k] = '{0}+{1}_{2}'.format(right_name[k], move , right_move[k])
                        else:
                            if right_snp_name[k]:
                                right_snp_name[k] = '{0}+{1}_{2}'.format(right_snp_name[k], move, alleles[left][0:1])
                            else:
                                right_snp_name[k] = '{0}_{1}'.format(move, alleles[left][0:1])
    for i in range(len(left_snp_name)):
        left_sequence[i] = ''.join(left_sequence[i])
        if left_snp_name[i]:
            left_name[i] = left_name[i] + ':' + left_snp_name[i]
    for i in range(len(right_snp_name)):
        right_sequence[i] = ''.join(right_sequence[i]) 
        if right_snp_name[i]:
            right_name[i] = right_name[i] + ':' + right_snp_name[i]
    return left_name, right_name, left_sequence, right_sequence



prefix = sys.argv[1]
offspring = sys.argv[2]
sample = sys.argv[3]
public_reference = sys.argv[4]
output_path = '{0}/{1}_{2}.fa'.format(prefix,offspring,sample)
name_output_path = '{0}/{1}_{2}.name.txt'.format(prefix,offspring,sample)            
reference = pysam.FastaFile(filename=public_reference)
variants = vcf.Reader(filename='{0}/{1}_{2}.hc.VQSR.wpost.denovo.phasing_plus.vcf.gz'.format(prefix,offspring,sample))
vcf_reader = vcf.Reader(filename='{0}/{1}_{2}.hc.VQSR.wpost.denovo.phasing_plus.vcf.gz'.format(prefix,offspring,sample))

#offspring = 'P6-1'
#sample = 'P6-female'
#output_path = 'bug.fa'
#name_output_path = 'bug.name.txt'
#reference = pysam.FastaFile(filename='../hg19.fa')
#variants = vcf.Reader(filename='test3.vcf.gz')
#vcf_reader = vcf.Reader(filename='test3.vcf.gz')
clear_file(output_path)
clear_file(name_output_path)

#chromosome = 'chr10'
#start = 105343462

chromosome = 'chr1'
start = 0
last_heter_pos = 0
heter_state = 0
sgRNA_align_length = 30
haplotype = {}
haplotype_name = {}
last_end = 0
haplotype_index = 0
for record in vcf_reader:
    record_chromosome = record.CHROM
    record_start = record.start
    record_end = record.end
#    print(record.samples)
    call = record.genotype(sample) 
    if call.is_het:
        if chromosome != record_chromosome:
#            print('next chromosome')
            print('{2}_{3}:{0}.{1}'.format(chromosome,record_chromosome,offspring,sample))
            if heter_state != 5:
                (left_name, right_name, left_sequence, right_sequence) = construct_haplotype(reference, variants, sample, chromosome, start)
                haplotype[haplotype_index+1] = left_sequence
                haplotype[haplotype_index+2] = right_sequence
                haplotype_name[haplotype_index+1] = left_name
                haplotype_name[haplotype_index+2] = right_name
                haplotype_index += 2
            else:
                (left_name_list, right_name_list, left_sequence_list, right_sequence_list) = construct_haplotype_state5(reference, variants, sample, chromosome, start)
                for i in range(len(left_name_list)):
                    haplotype[haplotype_index+i+1] = left_sequence_list[i]
                    haplotype_name[haplotype_index+i+1] = left_name_list[i]
                haplotype_index += len(left_sequence_list)
                for i in range(len(right_name_list)):
                    haplotype[haplotype_index+i+1] = right_sequence_list[i]
                    haplotype_name[haplotype_index+i+1] = right_name_list[i]
                haplotype_index += len(right_sequence_list)

            write_haplotype(output_path, haplotype)
            write_haplotype_name(name_output_path,haplotype_name)
            haplotype = {}
            haplotype_name = {}
            chromosome = record_chromosome
            start = 0
            last_heter_pos = 0
            heter_state = 0          
            last_end = 0      
        if call.phased:#or 'PGT' in str(call):
            if heter_state == 0:
                heter_state = 1
            elif heter_state == 1:
                heter_state = 1
            elif heter_state == 2:
                if record_start-last_heter_pos >= sgRNA_align_length and record_start > last_end :
                    end = record_start 
                    # left close and right open [)
                    (left_name, right_name, left_sequence, right_sequence) = construct_haplotype(reference, variants, sample, chromosome, start, end)
                    haplotype[haplotype_index+1] = left_sequence
                    haplotype[haplotype_index+2] = right_sequence
                    haplotype_name[haplotype_index+1] = left_name
                    haplotype_name[haplotype_index+2] = right_name
                    haplotype_index += 2

                    heter_state = 1
                    start = last_heter_pos
                else:
                    heter_state = 5                                                 
            elif heter_state == 5:
                if record_start-last_heter_pos >= sgRNA_align_length and record_start > last_end :
                    end = record_start 
                    # left close and right open [)
                    (left_name_list, right_name_list, left_sequence_list, right_sequence_list) = construct_haplotype_state5(reference, variants, sample, chromosome, start, end)
                    print("start" + str(start) + "\t" + "end" + str(end)) 
                    for i in range(len(left_name_list)):
                        haplotype[haplotype_index+i+1] = left_sequence_list[i]
                        haplotype_name[haplotype_index+i+1] = left_name_list[i]
                    haplotype_index += len(left_sequence_list)
                    for i in range(len(right_name_list)):
                        haplotype[haplotype_index+i+1] = right_sequence_list[i]
                        haplotype_name[haplotype_index+i+1] = right_name_list[i]
                    haplotype_index += len(right_sequence_list)
                        
                    heter_state = 1
                    start = last_heter_pos
                else:
                    heter_state = 5                       
            last_heter_pos = record_end
        else:
            if heter_state == 0:
                heter_state = 2 
            elif heter_state == 1:
                if record_start-last_heter_pos >= sgRNA_align_length and record_start > last_end :
                    end = record_start 
                    # left close and right open [)
                    (left_name, right_name, left_sequence, right_sequence) = construct_haplotype(reference, variants, sample, chromosome, start, end)
                    haplotype[haplotype_index+1] = left_sequence
                    haplotype[haplotype_index+2] = right_sequence
                    haplotype_name[haplotype_index+1] = left_name
                    haplotype_name[haplotype_index+2] = right_name
                    haplotype_index += 2
                    heter_state = 2
                    start = last_heter_pos
                else:
                    heter_state = 5
            elif heter_state == 2:
                if record_start-last_heter_pos >= sgRNA_align_length and record_start > last_end :
                    end = record_start 
                    # left close and right open [)
                    (left_name, right_name, left_sequence, right_sequence) = construct_haplotype(reference, variants, sample, chromosome, start, end)
                    haplotype[haplotype_index+1] = left_sequence
                    haplotype[haplotype_index+2] = right_sequence
                    haplotype_name[haplotype_index+1] = left_name
                    haplotype_name[haplotype_index+2] = right_name
                    haplotype_index += 2
                    heter_state = 2
                    start = last_heter_pos
                else:
                    heter_state = 5 
            elif heter_state == 5:
                if record_start-last_heter_pos >= sgRNA_align_length and record_start > last_end :
                    end = record_start 
                    # left close and right open [)
                    (left_name_list, right_name_list, left_sequence_list, right_sequence_list) = construct_haplotype_state5(reference, variants, sample, chromosome, start, end)
                    print("start" + str(start) + "\t" + "end" + str(end)) 
                    for i in range(len(left_name_list)):
                        haplotype[haplotype_index+i+1] = left_sequence_list[i]
                        haplotype_name[haplotype_index+i+1] = left_name_list[i]
                    haplotype_index += len(left_sequence_list)
                    for i in range(len(right_name_list)):
                        haplotype[haplotype_index+i+1] = right_sequence_list[i]
                        haplotype_name[haplotype_index+i+1] = right_name_list[i]
                    haplotype_index += len(right_sequence_list)
                    heter_state = 2
                    start = last_heter_pos
                else:
                    heter_state = 5
            last_heter_pos = record_end
    last_end = max(record_end,last_end)
else:
    print('end chromosome')
    print(chromosome)
    if heter_state != 5:
        (left_name, right_name, left_sequence, right_sequence) = construct_haplotype(reference, variants, sample, chromosome, start)
        haplotype[haplotype_index+1] = left_sequence
        haplotype[haplotype_index+2] = right_sequence
        haplotype_name[haplotype_index+1] = left_name
        haplotype_name[haplotype_index+2] = right_name
        haplotype_index += 2
    else:
        (left_name_list, right_name_list, left_sequence_list, right_sequence_list) = construct_haplotype_state5(reference, variants, sample, chromosome, start)
        for i in range(len(left_name_list)):
            haplotype[haplotype_index+i+1] = left_sequence_list[i]
            haplotype_name[haplotype_index+i+1] = left_name_list[i]
        haplotype_index += len(left_sequence_list)
        for i in range(len(right_name_list)):
            haplotype[haplotype_index+i+1] = right_sequence_list[i]
            haplotype_name[haplotype_index+i+1] = right_name_list[i]
        haplotype_index += len(right_sequence_list)

    write_haplotype(output_path, haplotype)       
    write_haplotype_name(name_output_path,haplotype_name)
    haplotype = {}
    haplotype_name = {}

