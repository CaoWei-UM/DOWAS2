#!/bin/bash

configure_file=${1}
check_dir()
{
	if [ ! -d "ref" ]; then
    	# Control will enter here if $DIRECTORY do not exists.
    		mkdir ref
	fi
}

check_varients_reference()
{
	if  !([ -f ${1} ] && [ -f ${1}.idx ]) ; then
		ref_name=$(basename ${1}) 
		if !([ -f ./ref/$ref_name ] && [ -f ./ref/$ref_name.idx ]); then
				cp ${1} ./ref/$ref_name
				$java -jar $gatk4 IndexFeatureFile -F $ref_file
		fi
		echo "$ref_file=./ref/$ref_name" >> $configure_file
	fi
}

check_public_sequence_reference()
{
	ref_prefix_name=$(basename $reference .fa)
	if !([ -f $reference ] && [ -f $reference.amb ] && [ -f  ${reference%.fa}.dict ] && [ -f $reference.ann ] && [ -f $reference.bwt ] && [ -f $reference.fai ] && [ -f $reference.pac ] && [ -f $reference.sa ] ) ; then
		ref_name=$(basename $reference)
		if !([ -f ./ref/$reference ] || [ -f ./ref/$reference.fa ]); then
				cp $reference ./ref/$ref_name
				$samtools faidx ./ref/$ref_name
				$bwa index ./ref/$ref_name
				$java -Xmx20g -jar $gatk4 CreateSequenceDictionary -R ./ref/$ref_name -O ./ref/$ref_prefix_name.dict
		fi
		echo "reference=./ref/$ref_prefix_name" >> $configure_file
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
check_dir
check_varients_reference $dbsnp
check_varients_reference $hapmap
check_varients_reference $OneKG_snp
check_varients_reference $OneKG_omni
check_varients_reference $OneKG_indel
check_varients_reference $Mills_indel
check_public_sequence_reference
check_sgRNA $guide_RNA_seq