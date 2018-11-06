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
parent=${3}
configure_file=${4}
output_dir=./${family_name}/output
metadata_dir=./${family_name}/metadata
log_dir=./${family_name}/log
source $configure_file
sgRNA_query=./ref/${guide_RNA_seq}_sgRNA_query.fa
################ function region ###############


select_parent_Variants(){
	$java -Xmx20g -jar ${gatk4} \
	     SelectVariants \
	     --exclude-filtered true  \
	     --sample-name ${parent} \
	     --remove-unused-alternates true \
	     --exclude-non-variants true \
	     -V ${metadata_dir}/mother_father_${offspring}.hc.VQSR.wpost.denovo.phased.vcf \
	     -O ${metadata_dir}/${offspring}_${parent}.hc.VQSR.wpost.denovo.vcf \
	     1> ${log_dir}/${offspring}.select_${parent}.log 2>&1 
}

phasing_plus(){

	    python ./phasing_plus.py ${metadata_dir} ${offspring} ${parent} 

}

index_Variants(){

	bgzip -c ${metadata_dir}/${offspring}_${parent}.hc.VQSR.wpost.denovo.phasing_plus.vcf > ${metadata_dir}/${offspring}_${parent}.hc.VQSR.wpost.denovo.phasing_plus.vcf.gz && \
	    tabix -p vcf ${metadata_dir}/${offspring}_${parent}.hc.VQSR.wpost.denovo.phasing_plus.vcf.gz 

}

construct_haplotypes(){

	    python ./construct_haplotypes_3.py ${metadata_dir} ${offspring} ${parent} ${reference}

}

mapping_to_hap(){

	    line_num_return=`wc -l ${metadata_dir}/${offspring}_${parent}.fa`
	    line_num=`echo ${line_num_return} | cut -d " " -f 1`    
	    line_split=`expr ${line_num} / 500 "*" 26 + 2` #cannot be /250*13, because it should be *(2N) to have even number lines in each file
	    split -l ${line_split} -d ${metadata_dir}/${offspring}_${parent}.fa ${metadata_dir}/${offspring}_${parent}.fa.part --verbose

	    for i in  00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
	    do
		mv ${metadata_dir}/${offspring}_${parent}.fa.part${i} ${metadata_dir}/${offspring}_${parent}.part${i}.fa && \
		   $bwa index ${metadata_dir}/${offspring}_${parent}.part${i}.fa && \
		    $bwa mem -T 23 -a ${metadata_dir}/${offspring}_${parent}.part${i}.fa ${sgRNA_query} > ${metadata_dir}/${offspring}_${parent}.part${i}.sam \
			2>${log_dir}/${offspring}_${parent}.part${i}_bwa.log  && \
		    $samtools view -bS -F 4 ${metadata_dir}/${offspring}_${parent}.part${i}.sam >${metadata_dir}/${offspring}_${parent}.part${i}.bam && \
		    $samtools sort ${metadata_dir}/${offspring}_${parent}.part${i}.bam -o ${metadata_dir}/${offspring}_${parent}.part${i}.sorted.bam &
		
	    done
		wait
}

hap_sam_merge(){

	    parted_bam=`ls ${metadata_dir}/${offspring}_${parent}.part*.sorted.bam`
	    $samtools merge -f -p -l 5 ${metadata_dir}/${offspring}_${parent}.sorted.bam ${parted_bam} 

}
hap_sam_to_ref_bed(){

	    python ./hap_sam_to_ref_bed.py ${metadata_dir} ${offspring} ${parent} 

}





    select_parent_Variants
    phasing_plus
    index_Variants
    construct_haplotypes
    mapping_to_hap
    hap_sam_merge
    hap_sam_to_ref_bed

