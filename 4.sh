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
chromosome=${2}
sample_gvcf_files=${3}
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



merge_samples_gvcfs_per_chromosome(){
	#merge samples from gvcf files per chromosome
	 $java -Xmx20g -jar $gatk4 CombineGVCFs -R "${reference}" ${sample_gvcf_files} -O "${metadata_dir}/${family_name}.hc.${chromosome}.g.vcf" 1> "${log_dir}/${family_name}.${chromosome}.combineGVCFs.log" 2>&1 
}

genotyping_gvcfs_per_chromosome(){
	$java -Xmx5g -jar ${gatk4} GenotypeGVCFs \
	 -R "${reference}" \
	 -V "${metadata_dir}/${family_name}.hc.${chromosome}.g.vcf" \
	 -O "${metadata_dir}/${family_name}.hc.${chromosome}.vcf" \
	 1>> "${log_dir}/${family_name}.${chromosome}.genotypeGVCFs.log" 2>&1 
}


################ main program ##################
source $configure_file
merge_samples_gvcfs_per_chromosome
genotyping_gvcfs_per_chromosome

