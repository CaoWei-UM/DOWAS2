# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 10:32:09 2018

@author: caowei
"""
import os
import sys
import time
import xml.etree.ElementTree as ET

chrom=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')
project_name=sys.argv[1]
if not os.path.exists('./{}/ref'.format(project_name)):
    os.makedirs('./{}/ref'.format(project_name))
tree = ET.parse('{}.xml'.format(project_name))
root = tree.getroot()
config='./{}/ref/config.txt'.format(project_name)
config_file=open(config, mode='w')
for depends in root.find('refdata_list'):
    config_file.writelines('{}=\'{}\'\n'.format(depends.tag,depends.text))
config_file.close()
member_list=[]
for seq_data in root.find('rawdata_list'):
    if(seq_data.tag not in locals().keys()):
        exec('{}=[]'.format(seq_data.tag))
        member_list.append(seq_data.tag)
    exec('{}.append({})'.format(seq_data.tag,seq_data.attrib))
if(os.system('./0.sh {}'.format(config))):
    print('FATAL ERROR: bulid reference error, please check your config xml file!')
    sys.exit(1)#建各种ref文件     
def make_ped_file(project_name,member_list):  
    if not os.path.exists('./{}/metadata'.format(project_name)):
        os.makedirs('./{}/metadata'.format(project_name))
    family_ped_file=open('./{}/metadata/{}.ped'.format(project_name,project_name), mode='w')	
    family_ped_file.writelines("5\tmother\t0\t0\t2\tfemale\n")
    family_ped_file.writelines("5\tfather\t0\t0\t1\tmale\n")
    for member in member_list: 
        if member !='father' and member != 'mother':
            offspring_ped_file=open('./{}/metadata/{}.ped'.format(project_name,member), mode='w')
            offspring_ped_file.writelines("5\tmother\t0\t0\t2\tfemale\n")
            offspring_ped_file.writelines("5\tfather\t0\t0\t1\tmale\n")
            family_ped_file.writelines("5\t{}\tfather\tfather\t5\tchild\n".format(member))
            offspring_ped_file.writelines("5\t{}\tfather\tfather\t5\tchild\n".format(member))
            offspring_ped_file.close()
    family_ped_file.close()
    return

make_ped_file(project_name,member_list)#还有ped文件
for member in member_list: 
    pid = os.fork()
    if pid == 0: 
        for lane in range(len(eval(member))): 
            os.system('./1.sh {} {} {} {} {} {}'.format(project_name,member,lane,eval(member)[lane]['R1'],eval(member)[lane]['R2'],config)) 
        os.system('./2.sh {} {} {}'.format(project_name,member,config))
        exit()
time.sleep(1)
os.wait()    

merge_vcf_list=[]
for chromosome in chrom:
    merge_vcf_list.append(" -I ./{}/metadata/{}.hc.{}.vcf".format(project_name,project_name,chromosome,config))
    pid = os.fork()
    if pid == 0: 
        gvcf_sample_list=[]
        for member in member_list: 
            os.system('./3.sh {} {} {} {}'.format(project_name,member,chromosome,config))#'建一个样本的一个chromosome的gvcf'
            gvcf_sample_list.append(" -V ./{}/metadata/{}_hc.{}.g.vcf".format(project_name,member,chromosome,config)) 
        sample_gvcfs=''.join(gvcf_sample_list)
        os.system('./4.sh {} {} \'{}\' {}'.format(project_name,chromosome,sample_gvcfs,config))#对每个chromosome混合样本并合并每个chromosome的vcf
        exit()

time.sleep(1)
os.wait()  
merge_vcfs=''.join(merge_vcf_list)        
os.system('./5.sh {} \'{}\' {}'.format(project_name,merge_vcfs,config))#处理混合成的一个vcf

for member in member_list: 
    pid = os.fork()
    if pid == 0:
        if member !='father' and member != 'mother':
            os.system('./6.sh {} {} {}'.format(project_name,member,config))#处理每个子代的数据
            for parent in ('mother','father'):
                os.system('./7.sh {} {} {} {}'.format(project_name,member,parent,config))#处理每个子代的数据与单个亲本的数据
            os.system('./8.sh {} {} {}'.format(project_name,member,config))#最后对比varients和germline所得的区间,得到最终potentional offtarget
        exit()
time.sleep(1)
os.wait()        
    
    
    