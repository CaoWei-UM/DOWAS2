# DOWAS2 manual
### update date: 2018/11/6
### version 0.1.2

## Introduction
DOWAS2 is a software made of python script and linux shell script. It is made for find potential off-target sites from reference and germline sequencing data in human. It consists of 5 pl scripts and 9 shell script files, run.py is the main and enter script. All information you need to input should be written in a xml file like example.xml and run ‘python run.py example’ in your computer. in that case, DOWAS2 will open ‘/your/current/path/example.xml’ and generate a folder in current path with same name with xml file (in that case it names ‘example’) and put all results in this folder.
You need to install some software and download some reference files. You should write them in the xml file before using DOWAS2, those softwares and reference with recommend version are list at the end of this manual.
Xml file contain a ‘rawdata_list’ subtag, it at least contains a father tag and a mother tag which contain raw sequencing data from parents in germline.
Final result will be saved in /your/path/example/metadata/${offspring}.hc.VQSR.wpost.denovo.phased.indel.casEdited.vcf, you will get multiple offspring results if you have more than one offspring.
## The advantages of DOWAS2 compare with DOWAS1
#### Shortages of DOWAS1
1. DOWAS1 use public reference(e.g. hg19) as mapping reference, that may lead to illusory off-target sites due to private variants in individuals.
2. DOWAS1 use a hard filter (mismatch<= 6 with PAM=NGG and mismatch<= 4 with PAM=NAG) to judge if a variant is an off-target site.
3. DOWAS1 cannot detect large deletions (deletions larger than 150bp in 20 kilo bases around on-target site). 
4. DOWAS1 doesn’t report SNPs due to high false positive ratio.
5. DOWAS1 cannot demonstrate the coverage at each potential off-target region directly.
6. DOWAS1 doesn't consider hotspots predicted by experimental methods.
7. The scalability of DOWAS1 is not good.
#### DOWAS2
DOWAS2 have a better software stracture and can overcome shortage 1,3,5,6,7 in DOWAS1.
## Getting started
```
cd /your/path
tar czvf DOWAS2v0.1.2.tar.gz 
cd DOWAS2v0.1.2
# Write your own config file, you can refer to example.xml
python run.py example
```
## Availability
DOWAS2 is released under MIT license. The latest version of DOWAS2 can be downloaded at www.sustc-genome.edu.cn. and github.



## Software and database dependencies
fastp 0.19.4

samtools-1.9

jre-1.8.0_181

GATK-3.8

GATK-4.0.4.0

bwa-0.7.17

vcftools-0.1.15

bedtools-v2.27.1

mosdepth-0.2.3

dbsnp database

OneKG_indel database

Mills_indel database

Hapmap database

OneKG_omni database

OneKG_snp database

Human genome sequence reference
## Seeking help
If you have questions about DOWAS2, you may send the questions to chenyr@mail.sustc.edu.cn or chenkj@mail.sustc.edu.cn or 11749245@mail.sustc.edu.cn or caow@mail@sustc.edu.cn . You may also ask questions in forums such as BioStar and SEQanswers.

## Update info
#### DOWAS2v0.1.2
fix some bug, add reference and software check before run main program.

## Author
DOWAS2 is written by Chen Kaijing, Chen Yangran, ChenRui and Cao Wei.
## Reference
•	Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] 

•	Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

•	Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

•	Poplin R, Ruano-Rubio V, DePristo M A, et al. Scaling accurate genetic variant discovery to tens of thousands of samples[J]. bioRxiv, 2017: 201178.
