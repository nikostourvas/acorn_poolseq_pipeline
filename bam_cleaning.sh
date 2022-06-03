#!/bin/bash

BASE="/home/tourvasn/ngs_training"

# reference genome
REF=$BASE/data/reference/Qrob_PM1N.fa

# declare variables
IND=$1
FORWARD=$BASE/results/trimmed_fastq/${IND}_R_1.trim.fastq.gz
REVERSE=$BASE/results/trimmed_fastq/${IND}_R_2.trim.fastq.gz
OUTPUT=$BASE/results/align/${IND}

samtools view -@ 15 -bh -q 20 ${OUTPUT}.sort.bam > ${OUTPUT}.sort.Q20.bam
samtools rmdup ${OUTPUT}.sort.Q20.bam ${OUTPUT}.sort.Q20.nodup.bam
samtools index -@ 15 ${OUTPUT}.sort.Q20.nodup.bam
samtools flagstat ${OUTPUT}.sort.Q20.nodup.bam > ${OUTPUT}.sort.Q20.nodup.stats
rm ${OUTPUT}${IND}.raw.bam ${OUTPUT}${IND}.sort.bam ${OUTPUT}${IND}.sort.Q20.bam
