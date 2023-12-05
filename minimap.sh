#!/bin/bash

#Lars Littmann
#2023/12/05
#Mapping script for paired-end short reads using minimap2 as aligner. 

#declare variables
#Directories mentioned below need to be updated manually. There is no practical way to take this all as command line arguments.
SAMPLE=${1}
REF=../TEST_LL/REFERENCE/Qrob_PM1N_Organelles.fa
FORWARD=../mapping_tests/${IND}_1.trim.fastq.gz
REVERSE=../mapping_tests/${IND}_2.trim.fastq.gz
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}"
OUTPUT=../mapping_tests/results/${IND}_minimap
THREADS=15
MEMORY=$((${THREADS}*4))

#minimap2 arguments:
# -ax sr align short reads
# -t number of threads
# -I available memory. Determined at 4GB per thread. 

#samtools view arguments:
#-h include header
#-b output a bam file

minimap2 -ax sr -R ${RG} -t ${THREADS} -I ${MEMORY}G ${REF} ${FORWARD} ${REVERSE} 2> ${OUTPUT}.mimimap2.err \
| samtools view --threads ${THREADS} -h -b -o ${OUTPUT}.raw.bam 2> ${OUTPUT}.sam.view.err

#gather statistics
samtools flagstat -@ 1 ${OUTPUT}.raw.bam > ${OUTPUT}.raw.flagstat
