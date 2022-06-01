#!/bin/bash

DATA="/home/tourvasn/ngs_training/results/align/Plomion"
RESULTS="/home/tourvasn/ngs_training/results/VCF"
REF="/home/tourvasn/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa"
	
/programs/samtools/bin/samtools mpileup -B -q 1 -f $REF $DATA/*_sort.rmd.bam > $RESULTS/cohort.mpileup
java -jar /programs/VarScan.jar mpileup2snp $RESULTS/cohort.mpileup \
	--vcf-sample-list inds \
	--min-coverage 10 --min-var-freq 0.0125 --p-value 0.05 \
	--output-vcf 1 > $RESULTS/cohort.varScan.snp.vcf
