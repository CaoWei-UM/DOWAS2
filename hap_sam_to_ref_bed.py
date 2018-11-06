import pysam
import pandas as pd
from collections import defaultdict
import sys

def get_bed_from_bam(sgRNA_hap_bam,bedfile,hap_name_dict):
    for read in sgRNA_hap_bam:
        if read.is_reverse :
            strand = '-'
        else:
            strand = '+'
        name = str(read.reference_name)
        position = int(read.reference_start)
#        query_sequence = str(read.query_sequence)
        (mismatch,pam,_) = str(read.query_name).split('_')
        (chromosome,start,end) = bed_reference_to_haplotype(name, position, position + 23,hap_name_dict)
        if (chromosome == "chr10" and start == 105356863 ):
            print(name)
        bedfile['chromosome'].append(chromosome)
        bedfile['start'].append(start)
        bedfile['end'].append(end)
#        bedfile['sequence'].append(query_sequence)
        bedfile['mismatch'].append(mismatch)
        bedfile['strand'].append(strand)
        bedfile['PAM'].append(pam)
    return bedfile

def hap_name_to_dict(hap_name_site):
    hap_name_dict={}
    f = open(hap_name_site,"r")
    lines = f.readlines()
    f.close()
    for line in lines:
        (index,hap_name)=line.split('>>')
        hap_name_dict[index]=hap_name
    return hap_name_dict

def bed_reference_to_haplotype(name, haplotype_start, haplotype_end,hap_name_dict):
    hap_name = hap_name_dict[name]
    name_split_list = [str(i) for i in hap_name.split(':')]
#    print(name_split_list)
    haplotype_chromosome = name_split_list[0]
    reference_start = int(name_split_list[1])
    indel_record_s = name_split_list[2]
    indel_sum = 0
    start_state = 0
    end_state = 0
    last_indel_length = 0
    last_haplotype_move = 0
    last_move = 0
    indel_name = indel_record_s.split('+')
    for indel_record in indel_name:
        (move, indel_length) = [int(i) for i in indel_record.split('_') ]
        haplotype_move = move - indel_sum     
        if start_state == 0:
            if haplotype_start <= haplotype_move:                                    
                if haplotype_start <= last_haplotype_move - last_indel_length and last_indel_length < 0:
                    reference_start_move = last_move + 1 
                else:
                    reference_start_move = haplotype_start + indel_sum 
                start_state = 1
        if end_state == 0:
            if haplotype_end <= haplotype_move:                                    
                if haplotype_end <= last_haplotype_move - last_indel_length and last_indel_length < 0:
                    reference_end_move = last_move + 1
                else:
                    reference_end_move = haplotype_end + indel_sum 
                end_state = 1
        if start_state == 1 and end_state == 1:
            break
        last_move = move
        last_haplotype_move = haplotype_move
        last_indel_length = indel_length
        indel_sum += indel_length
        if start_state == 0:
            if haplotype_start <= last_haplotype_move - last_indel_length and last_indel_length < 0:
                reference_start_move = last_move + 1
            else:
                reference_start_move = haplotype_start + indel_sum
        if end_state == 0:
            if haplotype_end <= last_haplotype_move - last_indel_length and last_indel_length < 0:
                reference_end_move = last_move + 1
            else:
                reference_end_move = haplotype_end + indel_sum
    return haplotype_chromosome, reference_start_move + reference_start, reference_end_move + reference_start

#sgRNA = 'CCR5'
#haplotype = 'P5-male.5'


prefix = sys.argv[1]
offspring = sys.argv[2]
sample = sys.argv[3]

bedfile_path = '{0}/{1}_{2}.bed'.format(prefix,offspring,sample)
sgRNA_hap_bam = pysam.AlignmentFile(filename='{0}/{1}_{2}.sorted.bam'.format(prefix,offspring,sample))
hap_name_site = '{0}/{1}_{2}.name.txt'.format(prefix,offspring,sample)

hap_name_dict = hap_name_to_dict(hap_name_site)
bedfile = defaultdict(list)
bedfile = get_bed_from_bam(sgRNA_hap_bam,bedfile,hap_name_dict)
bedfile = pd.DataFrame(bedfile)
bedfile.to_csv(bedfile_path, sep = '\t', index = False, header = False)


