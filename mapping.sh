#!/bin/bash

# create output directory
mkdir -p ../results/align

# reference genome
REF=../data/reference/Qrob_PM1N.fa

# declare variables
IND=${1}
FORWARD=../results/trimmed_fastq/${IND}_1.trim.fastq.gz
REVERSE=../results/trimmed_fastq/${IND}_2.trim.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=../results/align/${IND}

# then align, filter for quality, remove duplicates, and sort
# -t how many cores to use PER SAMPLE
bwa mem -R ${RG} -M -t 16 ${REF} ${FORWARD} ${REVERSE} > ${OUTPUT}.raw.sam

# -@ how many cores to use PER SAMPLE
samtools view -@ 16 -b ${OUTPUT}.raw.sam > ${OUTPUT}.raw.bam
rm ${OUTPUT}.raw.sam
# sort reads in the BAM according to their names so that pairs are placed one below the other
samtools sort -n -@ 16 ${OUTPUT}.raw.bam > ${OUTPUT}.name.sort.bam
rm ${OUTPUT}.raw.bam
# fix read-mates so that they both have the same sets of attributesfor the subsequent preprocessing
samtools fixmate -m ${OUTPUT}.name.sort.bam ${OUTPUT}.fixmate.sort.bam
# gather statistics
samtools flagstat ${OUTPUT}.fixmate.sort.bam > ${OUTPUT}.sort.flagstat