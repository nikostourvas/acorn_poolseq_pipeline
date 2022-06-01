#!/bin/bash

BASE="/home/tourvasn/ngs_training/results/align/Plomion"
IND=$1

java -jar /programs/picard.jar \
	MarkDuplicates REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	INPUT=$BASE/${IND}_sort.bam \
	OUTPUT=$BASE/${IND}_sort.rmd.bam \
	METRICS_FILE=$BASE/${IND}_sort.rmd.bam.metrics
 
# Index all bam files again
samtools index -@ 10 $BASE/*.rmd.bam
