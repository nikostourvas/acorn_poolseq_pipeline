#!/bin/bash

# create output directory
mkdir -p ../results/VCF

# declare variables
DATA=../results/align
RESULTS=../results/VCF
REF=../data/reference/Qrob_PM1N.fa
REGION=$1
	
samtools mpileup -B -q 1 -r $REGION -f $REF $DATA/*.sort.Q20.markdup.bam \
    > $RESULTS/$REGION.mpileup

java -jar /usr/share/java/varscan.jar mpileup2snp $RESULTS/$REGION.mpileup \
	--vcf-sample-list inds \
	--min-coverage 8 --min-var-freq 0.0055 --p-value 0.05 \
	--output-vcf 1 > $RESULTS/$REGION.varScan.snp.vcf

java -jar /usr/share/java/varscan.jar mpileup2indel $RESULTS/$REGION.mpileup \
	--vcf-sample-list inds \
	--min-coverage 8 --min-var-freq 0.0055 --p-value 0.1 \
	--output-vcf 1 > $RESULTS/$REGION.varScan.indel.vcf

#rm $RESULTS/$REGION.mpileup