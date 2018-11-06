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

bed_merge(){

	awk '{print $1"\t"$2"\t"$3}' ${metadata_dir}/${offspring}_mother.bed >${metadata_dir}/${offspring}_parents.bed && \
 	    awk '{print $1"\t"$2"\t"$3}' ${metadata_dir}/${offspring}_father.bed >>${metadata_dir}/${offspring}_parents.bed && \
	    cat ${metadata_dir}/${offspring}_parents.bed|sort|uniq >${metadata_dir}/${offspring}_parents.uniq.bed 

}
remove_regions_that_intersect_with_N_gaps(){

	$bedtools intersect -v -a ${metadata_dir}/${offspring}_parents.uniq.bed -b ./gaps_hg19.sorted.bed > ${metadata_dir}/${offspring}_parents.uniq.Nremoved.bed 

}
calculate_off_target_depth(){

	$mosdepth -n --by ${metadata_dir}/${offspring}_parents.uniq.Nremoved.bed --thresholds 1,2,9,12 ${metadata_dir}/${offspring}_parents.uniq.Nremoved ${metadata_dir}/${offspring}.sorted.markdup.BQSR.bam 

}
intersect_to_find_off_target(){

	    $bedtools intersect -wa -a ${metadata_dir}/${offspring}.hc.VQSR.wpost.denovo.phased.indel.vcf -b ${metadata_dir}/${offspring}_parents.uniq.Nremoved.bed > ${metadata_dir}/${offspring}.hc.VQSR.wpost.denovo.phased.indel.casEdited.vcf 

}

source $configure_file
    bed_merge
    remove_regions_that_intersect_with_N_gaps
    calculate_off_target_depth
    intersect_to_find_off_target

