#!/bin/bash

BASE="/home/tourvasn/ngs_training"

# create output directory
mkdir $BASE/results/align

# align a single individual
REF=$BASE/data/reference/Qrob_PM1N.fa

# declare variables
IND=$1
FORWARD=$BASE/results/trimmed_fastq/${IND}_R_1.trim.fastq.gz
REVERSE=$BASE/results/trimmed_fastq/${IND}_R_2.trim.fastq.gz
OUTPUT=$BASE/results/align/${IND}_sort.bam

# then align, filter for quality, remove duplicates, and sort
# -t how many cores to use PER SAMPLE
# -@ how many cores to use PER SAMPLE
bwa mem -M -t 15 $REF $FORWARD $REVERSE | \
samtools view -@ 15 -q 20 -b | \
samtools rmdup | \
samtools sort -@ 15 -T ${IND} > $OUTPUT

# statistics about sorted bam file
samtools flagstat $OUTPUT > $BASE/results/align/${IND}_sort.bam.flagstats
