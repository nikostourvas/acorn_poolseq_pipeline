#!/bin/bash

BASE="/home/tourvasn/ngs_training/results/align/Plomion/"

# -j number of jobs to run in parallel
cat inds | parallel --verbose -j 10 \
	picard \
	MarkDuplicates REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	INPUT=$BASE{}.bam \
	OUTPUT=$BASE{}.rmd.bam \
	METRICS_FILE=$BASE{}.rmd.bam.metrics
 
# Index all bam files again
samtools index $BASE*.rmd.bam
