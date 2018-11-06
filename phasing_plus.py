#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 01:35:10 2018

@author: chenyr
"""

import vcf
import sys
from collections import defaultdict



#@profile
def phasing():
    sample = sys.argv[3]
    offspring = sys.argv[2]
    prefix = sys.argv[1]
    vcf_reader = vcf.Reader(filename='{0}/{1}_{2}.hc.VQSR.wpost.denovo.vcf'.format(prefix,offspring,sample))
    #sample = "P6-male"
    #vcf_reader = vcf.Reader(filename="../metadata/P6-1/P6-1_P6-male.hc.VQSR.wpost.denovo.vcf")






    mapping_strand_dict = defaultdict(int)
    phasing_plus = []
    record_list = []
#    i = 0
#    record = next(vcf_reader)
    for i,record in enumerate(vcf_reader):
#    while record:
        record_list.append(record)
        call = record.genotype(sample)
        if 'PGT' in str(call):
            PID = str(call['PID'])
            (PGT_left, PGT_right) = [int(j) for j in call['PGT'].split("|")]
            if call.phased:
                (left, right) = [int(j) for j in call['GT'].split("|")]
                if mapping_strand_dict[PID] != 5: # mapping strand(0,1,2,5) 0:default 1:ALL GT=PGT 2:ALL GT=reverse(PGT) 5:some GT=PGT,some GT=reverse(PGT)
                    if (left, right) == (PGT_left, PGT_right):
                        if mapping_strand_dict[PID] == 2:
                            mapping_strand_dict[PID] == 5
                        elif mapping_strand_dict[PID] == 0:
                            mapping_strand_dict[PID] == 1
                    else:    
                        if mapping_strand_dict[PID] == 1:
                            mapping_strand_dict[PID] == 5
                        elif mapping_strand_dict[PID] == 0:
                            mapping_strand_dict[PID] == 1
            else:
                phasing_plus.append(i)
#        if i%1000 == 0:         
#            print('read_{0}'.format(i))
#        i += 1
#        record = next(vcf_reader)


    to_phasing_set=set(phasing_plus)
#    vcf_reader = vcf.Reader(filename='{0}/{1}/{1}_{2}.hc.VQSR.wpost.denovo.vcf'.format(prefix,offspring,sample))
    vcf_writer = vcf.Writer(open('{0}/{1}_{2}.hc.VQSR.wpost.denovo.phasing_plus.vcf'.format(prefix,offspring,sample), 'w'), vcf_reader)
    #vcf_reader = vcf.Reader(filename="../metadata/P6-1/P6-1_P6-male.hc.VQSR.wpost.denovo.vcf")
    #vcf_writer = vcf.Writer(open('../metadata/P6-1/P6-1_P6-male.hc.VQSR.wpost.denovo.phasing_plus.vcf', 'w'), vcf_reader)
           
            
    for i,record in enumerate(record_list):
        call=record.genotype(sample)
        if i in to_phasing_set:
            PID = str(call['PID']) 
            if mapping_strand_dict[PID] == 2:
                record.genotype(sample).data = call.data._replace(GT=str(call['PGT'][::-1]))
            elif mapping_strand_dict[PID] == 1:
                record.genotype(sample).data = call.data._replace(GT=str(call['PGT']))
            elif mapping_strand_dict[PID] == 0 and 'PGT' in str(call):
                record.genotype(sample).data = call.data._replace(GT=str(call['PGT']))
        vcf_writer.write_record(record)
#        if i%1000 == 0:
#            print('write_{0}'.format(i))    
            
            
            
phasing()
