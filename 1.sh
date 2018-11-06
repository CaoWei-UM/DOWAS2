#!/bin/bash

######## input ##########
# $1 family name(project name)
# $2 sample name
# $3 lane
# $4 Path to pair1 sequence rawdata
# $5 Path to pair2 sequence rawdata
family_name=${1}
sample_name=${2}
lane=${3}
rawdata_r1=${4}
rawdata_r2=${5}
configure_file=${6}
output_dir=./${family_name}/output
metadata_dir=./${family_name}/metadata
log_dir=./${family_name}/log


##################################### constant #################################

#RG="@RG\tID:${sample_name}\tLB:${sample_name}L\tPL:illumina\tSM:${sample_name}\tPU:${sample_name}"
#reference=/1/public_resources/humanGenome/hg19/hg19_bwa/hg19_bwa
################################### function ####################################

check_dir()
{
	if [ ! -d "${1}" ]; then
    	# Control will enter here if $DIRECTORY do not exists.
    		mkdir "${1}"
	fi
}

quality_control()
{   
    # improve the sequence data quality
    $fastp -w 8 -c -h "${output_dir}/${sample_name}.${lane}.qc.html" -j "${output_dir}/${sample_name}.${lane}.qc.json" -i "${rawdata_r1}" -o "${metadata_dir}/${sample_name}.r1.qc.fq.gz" -I "${rawdata_r2}" -O "${metadata_dir}/${sample_name}.r2.qc.fq.gz"
}

mapping()
{
    # Align the paired reads to reference genome using bwa mem && samtobam with samtools view -bS
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>align to reference\n" 
    $bwa mem -t 20 -M -Y\
	-R "@RG\tID:${sample_name}\tPL:ILLUMINA\tPU:HISEQ\tLB:${sample_name}library\tSM:${sample_name}" \
	"${reference}" \
	"${metadata_dir}/${sample_name}.r1.qc.fq.gz" \
	"${metadata_dir}/${sample_name}.r2.qc.fq.gz" \
	2>> "${log_dir}/${sample_name}.bwa_mem.log" \
	| $samtools view -@ 20 -Sb - > "${metadata_dir}/${sample_name}.bam" \
		   2>> "${log_dir}/${sample_name}.${lane}.sam_to_bam.log"

}

sortbam(){
    # sort bam     
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>bam to sorted"
    $samtools sort -m 4G -@ 16 \
        -o "${metadata_dir}/${sample_name}.${lane}.sorted.bam" \
        "${metadata_dir}/${sample_name}.bam" \
        1>> "${log_dir}/${sample_name}.${lane}.sort.log" 2>&1

	echo "${metadata_dir}/${sample_name}.${lane}.sorted.bam" >> ${metadata_dir}/${sample_name}.merge_list.txt
}

########################################### main ########################################


check_dir "${output_dir}"
check_dir "${metadata_dir}"
check_dir "${log_dir}"
source $configure_file

quality_control
mapping
sortbam
rm ${metadata_dir}/${sample_name}.bam
rm ${metadata_dir}/${sample_name}.r1.qc.fq.gz
rm ${metadata_dir}/${sample_name}.r2.qc.fq.gz