#!/bin/bash

RESULTS=../results/VCF

java -jar /programs/VarScan.jar mpileup2indel $RESULTS/cohort.mpileup \
	--min-coverage 10 --min-var-freq 0.0100 --p-value 0.10 \
	> $RESULTS/cohort.varScan.indel
	
java -jar VarScan.jar filter $RESULTS/cohort.varScan.snp \
	--indel-file $RESULTS/cohort.varScan.indel \
	--output-file $RESUTLS/cohort.varScan.snp.filter

#java –jar VarScan.jar filter sample.varScan.indel -–min-reads2 4 –-min-var-freq 0.15 –-p-value 0.05 –-output-file sample.varScan.indel.filter
