#!/bin/bash

UNTRIMMED=/home/tourvasn/ngs_training/data/untrimmed_fastq/
TRIMMED=/home/tourvasn/ngs_training/results/trimmed_fastq/

for infile in ${UNTRIMMED}*1_001.fastq.gz
 do
   base=$(basename ${infile} 1_001.fastq.gz)
   java -jar /programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${infile} ${UNTRIMMED}${base}2_001.fastq.gz \
                ${TRIMMED}${base}_1.trim.fastq.gz ${TRIMMED}${base}_1un.trim.fastq.gz \
                ${TRIMMED}${base}_2.trim.fastq.gz ${TRIMMED}${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15 \
		-threads 20
 done

