#!/bin/bash

######## input ##########
# $1 sample name
# $2 Path to pair1 sequence rawdata
# $3 Path to pair2 sequence rawdata
# $4 Path to output directory, it will be created if it does not exist.
# $5 Path to metadata directory, it will be created if it does not exist.
# $6 Path to log directory, it will be created if it does not exist.
family_name=${1}
sample_name=${2}
configure_file=${3}
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

mergebam(){
    # merge bam     
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>bam to merge"
    $samtools merge -b  \
        "${metadata_dir}/${sample_name}.merge_list.txt" \
        "${metadata_dir}/${sample_name}.merged.bam" \
        1>> "${log_dir}/${sample_name}.merge.log" 2>&1
	rm ${metadata_dir}/${sample_name}.*.sorted.bam
		
}

replaceRG(){
    # ReplaceReadGroup
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Replace ReadGroup"
    $java -Xmx20g -jar ${picard}  AddOrReplaceReadGroups \
        I="${metadata_dir}/${sample_name}.merged.bam" \
        O="${metadata_dir}/${sample_name}.sorted.RG.bam" \
		RGID="${sample_name}" \
        RGLB="${sample_name}L" \
        RGPL=illimina \
        RGPU=run \
        RGSM="${sample_name}" \
        1>> "${log_dir}/${sample_name}.RG.log" 2>&1
}

markdup(){
    # Mark duplicate reads
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Mark duplicates" 
    $java -Xmx20g -jar ${picard} MarkDuplicates \
	   I="${metadata_dir}/${sample_name}.sorted.RG.bam" \
	   O="${metadata_dir}/${sample_name}.sorted.RG.markdup.bam" \
	   M="${metadata_dir}/${sample_name}_markdup.metrics" \
           CREATE_INDEX=true \
	   REMOVE_DUPLICATES=true \
	   1>> "${log_dir}/${sample_name}.markdup.log" 2>&1
}

BQSR()
{
    # BaseRecalibration Calculate
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BaseRecalibration_Calculate"
    $java -Xmx50g -jar ${gatk4} \
        BaseRecalibrator \
        -R "${reference}" \
        -I "${metadata_dir}/${sample_name}.sorted.RG.markdup.bam" \
        -O "${metadata_dir}/${sample_name}_BQSR.metrics" \
        --known-sites "${OneKG_indel}" \
        --known-sites "${Mills_indel}" \
        --known-sites "${dbsnp}" \
        1>> "${log_dir}/${sample_name}.BQSR1.log" 2>&1

    # Print Recalibrated Reads
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BaseRecalibration_Print"
    $java -Xmx50g -jar ${gatk4} \
        ApplyBQSR \
        -bqsr "${metadata_dir}/${sample_name}_BQSR.metrics" \
        -I "${metadata_dir}/${sample_name}.sorted.RG.markdup.bam" \
        -O "${metadata_dir}/${sample_name}.sorted.RG.markdup.BQSR.bam" \
        1> "${log_dir}/${sample_name}.BQSR2.log" 2>&1
}



########################################### main ########################################
source $configure_file

mergebam
replaceRG
markdup
BQSR

