#!/bin/bash

# create output directory
mkdir -p ../results/align_Batch1

# declare variables
IND=${1}
REF=..reference/Qrob_PM1N.fa
FORWARD=../results/trimmed_fastq_Batch1/${IND}_1.trim.fastq.gz
REVERSE=../results/trimmed_fastq_Batch1/${IND}_2.trim.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=../results/align_Batch1/${IND}

# Align to reference genome and export raw bam file

# BWA mem arguments
# -M mark shorter split hits as secondary
# -t how many cores to use PER SAMPLE for mapping

# Samtools view arguments
# -h Include header in output
# -b output a bam file

# Be mindful to account for as many cores have been assigned to bwa and samtools 
# (e.g. 4c BWA + 1c samtools = 5 cores in total)
# samtools flagstat is run afterwards (much sorter time) utilizing all 5 cores

bwa mem -R ${RG} -M -t 4 ${REF} ${FORWARD} ${REVERSE} 2> ${OUTPUT}.bwa-mem.err \
    | samtools view -h -b ${OUTPUT}.raw.bam

# gather statistics
samtools flagstat -@ 5 ${OUTPUT}.raw.bam > ${OUTPUT}.raw.flagstat