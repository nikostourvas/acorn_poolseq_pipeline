#!/bin/bash

RESULTS="/home/tourvasn/ngs_training/results/align/Plomion/"
REF="/home/tourvasn/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa"

# -j number of jobs to run in parallel
cat inds | parallel --verbose -j 10 \
	samtools mpileup -f $REF $RESULTS{}_sort.bam > $RESULTS.mpileup	
