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
offspring=${2}
configure_file=${3}
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

select_family_per_offspring_by_GATK(){

	$java -Xmx20g -jar ${gatk4} \
             SelectVariants \
	     --exclude-filtered true  \
	     --sample-name mother \
	     --sample-name father \
	     --sample-name ${offspring} \
	     --remove-unused-alternates true \
	     --exclude-non-variants true \
	     -V ${metadata_dir}/${family_name}.hc.VQSR.wpost.denovo.vcf \
	     -O ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.gatk.vcf \
	     1> "${log_dir}/${offspring}.select.log" 2>&1 

}

phase_by_transmission(){


	 $java -Xmx20g -jar ${gatk3_8} \
	      -T PhaseByTransmission \
	      -R ${reference} \
	      -V ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.gatk.vcf \
	      -ped ${metadata_dir}/${offspring}.ped \
	      -o ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.phased.vcf \
	      1> "${log_dir}/${offspring}.gatk_phase.log" 2>&1 

}

select_offspring_indel_by_GATK(){

	    grep '#' ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.phased.vcf >${metadata_dir}/${offspring}.hc.VQSR.wpost.denovo.phased.indel.vcf && \
		grep "hiConfDeNovo" ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.phased.vcf | \
		    grep "${offspring}" \
			 >> ${metadata_dir}/${offspring}.hc.VQSR.wpost.denovo.phased.indel.vcf 
}


source $configure_file
    select_family_per_offspring_by_GATK
    phase_by_transmission
    select_offspring_indel_by_GATK
    