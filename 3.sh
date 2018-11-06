#!/bin/bash
## the second part of dowas, which assuming that you have produced gvcf of each sample with dowas1.sh

######## parameter #########
######## input ########
# $1 sample name
# $2 samples join by female,mother,and offspring list
# $3 Path to output directory, it will be created if it does not exist.
# $4 Path to metadata directory, it will be created if it does not exist.
# $5 Path to log directory, it will be created if it does not exist.
family_name=${1}
sample_name=${2}
chromosome=${3}
configure_file=${4}
output_dir=./${family_name}/output
metadata_dir=./${family_name}/metadata
log_dir=./${family_name}/log

################ function region ###############

check_dir()
{
    if [ ! -d "${1}" ]; then
	# Control will enter here if $DIRECTORY do not exists.
	mkdir "${1}"
    fi
}

HC_GVCF_caller_per_chr()
{
    # Haplotype Caller GVCF mode per chromosome
	$java -Xmx5g -jar ${gatk4} \
	     HaplotypeCaller \
	     --emit-ref-confidence GVCF \
	     -R "${reference}" \
	     -I "${metadata_dir}/${sample_name}.sorted.RG.markdup.BQSR.bam" \
	     -L "${chromosome}" \
	     -O "${metadata_dir}/${sample_name}_hc.${chromosome}.g.vcf" \
	     1> "${log_dir}/${sample_name}.hc.${chromosome}.log" 2>&1  
}

################ main program ##################
source $configure_file
HC_GVCF_caller_per_chr
