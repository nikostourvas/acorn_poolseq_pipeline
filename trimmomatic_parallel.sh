#!/bin/bash

# create output directory
mkdir -p ../results/trimmed_fastq

# variables
IND=$1
UNTRIMMED=/mnt/data/untrimmed_fastq/
TRIMMED=../results/trimmed_fastq/
ADAPTER=/mnt/data/adapters/TruSeq3-PE-2.fa

java -jar /usr/share/java/trimmomatic-0.39.jar PE \
	${UNTRIMMED}${IND}_R1_001.fastq.gz ${UNTRIMMED}${IND}_R2_001.fastq.gz \
        ${TRIMMED}${IND}_1.trim.fastq.gz ${TRIMMED}${IND}_1un.trim.fastq.gz \
        ${TRIMMED}${IND}_2.trim.fastq.gz ${TRIMMED}${IND}_2un.trim.fastq.gz \
        ILLUMINACLIP:${ADAPTER}:2:40:15 \
        SLIDINGWINDOW:4:20 \
        MINLEN:25 \ 
	-threads 16 \
        2> ${TRIMMED}trim_out.log
