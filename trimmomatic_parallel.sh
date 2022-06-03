#!/bin/bash

# create output directory
mkdir ../results/trimmed_fastq

# variables
IND=$1
UNTRIMMED=/home/tourvasn/ngs_training/data/untrimmed_fastq/
TRIMMED=/home/tourvasn/ngs_training/results/trimmed_fastq/

java -jar /programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	${UNTRIMMED}${IND}_1.fastq.gz ${UNTRIMMED}${IND}_2.fastq.gz \
        ${TRIMMED}${IND}_1.trim.fastq.gz ${TRIMMED}${IND}_1un.trim.fastq.gz \
        ${TRIMMED}${IND}_2.trim.fastq.gz ${TRIMMED}${IND}_2un.trim.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15 \
	-threads 15
