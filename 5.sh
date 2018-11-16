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
merge_vcfs=${2}
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




merge_vcf(){
## merge all vcf in a family
	 $java -Xmx5g -jar ${gatk4} \
	     MergeVcfs \
	     ${merge_vcfs} \
	     -O ${metadata_dir}/${family_name}.hc.vcf \
	     1>> "${log_dir}/${family_name}.merge_vcfs.log" 2>&1
}




vqsr_snp_mode(){
## VQSR SNP mode
     $java -Xmx20g -jar ${gatk4} VariantRecalibrator \
	 -R "${reference}" \
	 -V ${metadata_dir}/${family_name}.hc.vcf \
	 -resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
	 -resource omini,known=false,training=true,truth=false,prior=12.0:${OneKG_omni} \
	 -resource 1000G,known=false,training=true,truth=false,prior=10.0:${OneKG_snp} \
	 -resource dbsnp,known=true,training=false,truth=false,prior=6.0:${dbsnp} \
	 -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	 -mode SNP \
	 --max-gaussians ${1} \
	 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	 --tranches-file ${metadata_dir}/${family_name}.hc.snps.tranches \
	 -rscript-file ${metadata_dir}/${family_name}.hc.snps.plots.R \
	 -O ${metadata_dir}/${family_name}.hc.snps.recal \
	 1>> "${log_dir}/${family_name}.VQSR_snp_1.log" 2>&1 && \
	 $java -Xmx5g -jar ${gatk4} ApplyVQSR \
	     -R "${reference}" \
	     -V ${metadata_dir}/${family_name}.hc.vcf \
	     --ts-filter-level 99.0 \
	     --tranches-file ${metadata_dir}/${family_name}.hc.snps.tranches \
	     -recal-file ${metadata_dir}/${family_name}.hc.snps.recal \
	     -mode SNP \
	     -O ${metadata_dir}/${family_name}.hc.snps.VQSR.vcf \
	     1>> "${log_dir}/${family_name}.VQSR_snp_2.log" 2>&1 && \
	echo "** SNPs VQSR done **"
}

vqsr_indel_mode(){
## VQSR INDEL mode
     $java -Xmx20g -jar ${gatk4} VariantRecalibrator \
	 -R "${reference}" \
	 -V ${metadata_dir}/${family_name}.hc.snps.VQSR.vcf \
	 -resource mills,known=true,training=true,truth=true,prior=12.0:"${Mills_indel}" \
	 -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	 -mode INDEL \
	 --max-gaussians ${1} \
	 -rscript-file ${metadata_dir}/${family_name}.hc.snps.indels.plots.R \
	 --tranches-file ${metadata_dir}/${family_name}.hc.snps.indels.tranches \
	 -O ${metadata_dir}/${family_name}.hc.snps.indels.recal \
	 1>> "${log_dir}/${family_name}.VQSR_indel_1.log" 2>&1 && \
	 $java -Xmx5g -jar ${gatk4} ApplyVQSR \
	     -R "${reference}" \
	     -V ${metadata_dir}/${family_name}.hc.snps.VQSR.vcf \
	     --ts-filter-level 99.0 \
	     --tranches-file ${metadata_dir}/${family_name}.hc.snps.indels.tranches \
	     -recal-file ${metadata_dir}/${family_name}.hc.snps.indels.recal \
	     -mode INDEL \
	     -O ${metadata_dir}/${family_name}.hc.VQSR.vcf \
	     1>> "${log_dir}/${family_name}.VQSR_indel_2.log" 2>&1 && \
	echo "** SNPs and Indels VQSR (${family_name}.hc.VQSR.vcf finish) done **"
}


calculate_genotype_posteriors(){
     $java -Xmx20g -jar ${gatk4} CalculateGenotypePosteriors \
	 -R ${reference} \
	 -V ${metadata_dir}/${family_name}.hc.VQSR.vcf \
	 -O ${metadata_dir}/${family_name}.hc.VQSR.wpost.vcf \
	 -ped ${metadata_dir}/${family_name}.ped \
	 1>> ${log_dir}/${family_name}.wpost.log 2>&1 && \
	echo "** CalculateGenotypePosteriors (${family_name}.hc.VQSR.wpost.vcf finish) done **"
}

variant_annotation(){
     $java -jar ${gatk3_8} \
	 -T VariantAnnotator \
	 -R ${reference} \
	 -V ${metadata_dir}/${family_name}.hc.VQSR.wpost.vcf \
	 -A PossibleDeNovo \
	 -ped ${metadata_dir}/${family_name}.ped \
	 -o ${metadata_dir}/${family_name}.hc.VQSR.wpost.denovo.vcf \
	 1>> "${log_dir}/${family_name}.denovo.log" 2>&1 && \
	echo "** PossibleDeNovo annotate (${family_name}.hc.VQSR.wpost.denovo.vcf finish) done **"
   
}


################ main program ##################
source $configure_file
    merge_vcf #done
	gaussians_threshold=6
	while true; do
		vqsr_snp_mode $gaussians_threshold
		[ ! -e ${metadata_dir}/${family_name}.hc.snps.VQSR.vcf ] || break
		gaussians_threshold=`expr $gaussians_threshold - 2`
	done
	gaussians_threshold=6
	while true; do
		vqsr_indel_mode $gaussians_threshold #done
		[ ! -e ${metadata_dir}/${family_name}.hc.VQSR.vcf ] || break
		gaussians_threshold=`expr $gaussians_threshold - 2`
	done
    calculate_genotype_posteriors #done
    variant_annotation #done
	