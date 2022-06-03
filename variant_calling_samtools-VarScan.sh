#!/bin/bash

DATA="/home/tourvasn/ngs_training/results/align"
RESULTS="/home/tourvasn/ngs_training/results/VCF"
REF="/home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa"

mkdir ../results/VCF
	
samtools mpileup -B -q 1 -f $REF $DATA/*.sort.Q20.nodup.bam > $RESULTS/cohort.mpileup
java -jar /programs/VarScan.jar mpileup2snp $RESULTS/cohort.mpileup \
	--vcf-sample-list inds \
	--min-coverage 8 --min-var-freq 0.0055 --p-value 0.05 \
	--output-vcf 1 > $RESULTS/Qrob_209-220.varScan.snp.vcf
