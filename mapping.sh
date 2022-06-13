#!/bin/bash

BASE="/home/tourvasn/ngs_training"

# create output directory
mkdir -p ${BASE}/results/align

# reference genome
REF=${BASE}/data/reference/Qrob_PM1N.fa

# declare variables
IND=${1}
FORWARD=${BASE}/results/trimmed_fastq/${IND}_1.trim.fastq.gz
REVERSE=${BASE}/results/trimmed_fastq/${IND}_2.trim.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=${BASE}/results/align/${IND}

# then align, filter for quality, remove duplicates, and sort
# -t how many cores to use PER SAMPLE
bwa mem -R ${RG} -M -t 15 ${REF} ${FORWARD} ${REVERSE} > ${OUTPUT}.raw.sam

# -@ how many cores to use PER SAMPLE
samtools view -@ 15 -b ${OUTPUT}.raw.sam > ${OUTPUT}.raw.bam
# sort reads in the BAM according to their names so that pairs are placed one below the other
samtools sort -n -@ 15 ${OUTPUT}.raw.bam > ${OUTPUT}.name.sort.bam
# fix read-mates so that they both have the same sets of attributesfor the subsequent preprocessing
samtools fixmate -m ${OUTPUT}.name.sort.bam ${OUTPUT}.fixmate.sort.bam
# gather statistics
samtools flagstat ${OUTPUT}.fixmate.sort.bam > ${OUTPUT}.sort.stats
