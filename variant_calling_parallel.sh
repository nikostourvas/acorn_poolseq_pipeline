#!/bin/bash

# create output directory
mkdir -p ../results/VCF

# declare variables
DATA="/home/tourvasn/ngs_training/results/align"
RESULTS="/home/tourvasn/ngs_training/results/VCF"
REF="/home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa"
REGION=$1
	
samtools mpileup -B -q 1 -r $REGION -f $REF $DATA/*.sort.Q20.markdup.bam \
    > $RESULTS/$REGION.mpileup

java -jar /programs/VarScan.jar mpileup2snp $RESULTS/$REGION.mpileup \
	--vcf-sample-list inds \
	--min-coverage 8 --min-var-freq 0.0055 --p-value 0.05 \
	--output-vcf 1 > $RESULTS/$REGION.varScan.snp.vcf

java -jar /programs/VarScan.jar mpileup2indel $RESULTS/$REGION.mpileup \
	--vcf-sample-list inds \
	--min-coverage 8 --min-var-freq 0.0055 --p-value 0.1 \
	--output-vcf 1 > $RESULTS/$REGION.varScan.indel.vcf

rm $RESULTS/$REGION.mpileup