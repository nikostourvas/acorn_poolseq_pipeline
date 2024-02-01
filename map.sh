#!/bin/bash

# create output directory
mkdir -p /data/genetics_tmp/results/mapped_reads

# declare variables
IND=${1}
REF=/mnt/reference/Qrob_PM1N_Organelles.fa
FORWARD=/data/genetics_tmp/results/fastp_dedup_trim/${IND}_1.trim.dedup.fastq.gz
REVERSE=/data/genetics_tmp/results/fastp_dedup_trim/${IND}_2.trim.dedup.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=/data/genetics_tmp/results/mapped_reads/${IND}
BWAMEM2=/usr/local/bin/bwa-mem2/bwa-mem2 #path to executable needs to be explicit

# Align to reference genome and export raw bam file

# BWA mem2 arguments
# -M mark shorter split hits as secondary
# -t how many cores to use PER SAMPLE for mapping

# Samtools view arguments
# -h Include header in output
# -b output a bam file

# Be mindful to account for as many cores have been assigned to bwa and samtools 
# (e.g. 2c BWA + 1c samtools = 3 cores in total)
# samtools flagstat is run afterwards 

${BWAMEM2} mem -R ${RG} -M -t 1 ${REF} ${FORWARD} ${REVERSE} 2> ${OUTPUT}.bwa-mem.err \
    | samtools view --threads 1 -h -b -o ${OUTPUT}.raw.bam 2> ${OUTPUT}.sam-view.err

# gather statistics
samtools flagstat -@ 1 ${OUTPUT}.raw.bam > ${OUTPUT}.raw.flagstat
