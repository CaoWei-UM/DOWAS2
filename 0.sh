#!/bin/bash

configure_file=${1}
check_dir()
{
	if [ ! -d "ref" ]; then
    	# Control will enter here if $DIRECTORY do not exists.
    		mkdir ref
	fi
}

check_software()
{
	if  [ ! -f ${1} ] ; then
		software_name=`cat $configure_file|grep =\'${1}\'$|cut -d = -f 1|head -1` #reflection get software name
		echo "software $software_name not found!">&2
		exit 1
	fi
}
check_varients_reference()
{
	ref_varients_name=`cat $configure_file|grep =\'${1}\'$|cut -d = -f 1|head -1` #reflection get ref_varients name
	if  [ -f  ${1} ]; then
		if  [ ! -f ${1}.idx ] ;then
			ref_name=$(basename ${1}) 
			if !([ -f ./ref/$ref_name ] && [ -f ./ref/$ref_name.idx ]); then
					cp ${1} ./ref/$ref_name
					$java -jar $gatk4 IndexFeatureFile -F ./ref/$ref_name
			fi
			echo "$ref_varients_name='./ref/$ref_name'" >> $configure_file
		fi
	else
		echo "varients reference $ref_varients_name not found!">&2
		exit 1
	fi
}

check_public_sequence_reference()
{
	if  [ -f  $reference ]; then
		ref_prefix_name=$(basename $reference .fa)
		if !( [ -f $reference.amb ] && [ -f  ${reference%.fa}.dict ] && [ -f $reference.ann ] && [ -f $reference.bwt ] && [ -f $reference.fai ] && [ -f $reference.pac ] && [ -f $reference.sa ] ) ; then
			if [ ! -f ./ref/$ref_prefix_name.fa ] ; then
					cp $reference ./ref/$ref_prefix_name.fa
					$samtools faidx ./ref/$ref_prefix_name.fa
					$bwa index ./ref/$ref_prefix_name.fa
					$java -Xmx20g -jar $gatk4 CreateSequenceDictionary -R ./ref/$ref_prefix_name.fa -O ./ref/$ref_prefix_name.dict
			fi
			echo "reference='./ref/$ref_prefix_name.fa'" >> $configure_file
		fi
	else
		echo "reference not found!">&2
		exit 1
	fi
}

check_sgRNA()
{
	if  [ ! -f ./ref/${1}_sgRNA_query.fa ]; then
		python ./create_sgRNA.py ${1} ./ref/${1}_sgRNA_query.fa
	fi
}

########################################### main ########################################
source $configure_file
check_software $gatk3_8
check_software $gatk4
check_software $picard
check_software $samtools
check_software $bwa
check_software $java
check_software $fastp
check_software $bedtools
check_software $mosdepth
check_dir
check_varients_reference $dbsnp
check_varients_reference $hapmap
check_varients_reference $OneKG_snp
check_varients_reference $OneKG_omni
check_varients_reference $OneKG_indel
check_varients_reference $Mills_indel
check_public_sequence_reference
check_sgRNA $guide_RNA_seq
