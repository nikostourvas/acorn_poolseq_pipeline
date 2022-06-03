#!/bin/bash

BASE="/home/tourvasn/ngs_training"

# create output directory
mkdir $BASE/results/align

# reference genome
REF=$BASE/data/reference/Qrob_PM1N.fa

# declare variables
IND=$1
FORWARD=$BASE/results/trimmed_fastq/${IND}_R_1.trim.fastq.gz
REVERSE=$BASE/results/trimmed_fastq/${IND}_R_2.trim.fastq.gz
OUTPUT=$BASE/results/align/${IND}

# then align, filter for quality, remove duplicates, and sort
# -t how many cores to use PER SAMPLE
# -@ how many cores to use PER SAMPLE
bwa mem -R "@RG\tID:${IND}\tPL:Illumina\tSM:${IND}" \
	-M -t 15 $REF $FORWARD $REVERSE | \
samtools view -@ 15 -b > ${OUTPUT}.raw.bam
samtools sort -@ 15 ${OUTPUT}.raw.bam -o ${OUTPUT}.sort.bam
samtools flagstat ${OUTPUT}.sort.bam > ${OUTPUT}.sort.stats
