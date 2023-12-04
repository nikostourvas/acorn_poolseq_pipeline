#!/bin/bash

# create output directory
mkdir -p ../results/align_Batch2

# declare variables
IND=${1}
REF=../TEST_LL/REFERENCE/Qrob_PM1N_Organelles.fa
FORWARD=../mapping_tests/${IND}_1.trim.fastq.gz
REVERSE=../mapping_tests/${IND}_2.trim.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=../mapping_tests/results/${IND}

# Align to reference genome and export raw bam file

# BWA mem arguments
# -M mark shorter split hits as secondary
# -t how many cores to use PER SAMPLE for mapping

# Samtools view arguments
# -h Include header in output
# -b output a bam file

# Be mindful to account for as many cores have been assigned to bwa and samtools 
# (e.g. 2c BWA + 1c samtools = 3 cores in total)
# samtools flagstat is run afterwards 

bwa mem -R ${RG} -M -t 1 ${REF} ${FORWARD} ${REVERSE} 2> ${OUTPUT}.bwa-mem.err \
    | samtools view --threads 15 -h -b -o ${OUTPUT}.raw.bam 2> ${OUTPUT}.sam-view.err

# gather statistics
samtools flagstat -@ 1 ${OUTPUT}.raw.bam > ${OUTPUT}.raw.flagstat
