#!/bin/bash

# declare variables
IND=${1}
OUTPUT=../results/align/${IND}

# Produce filtered BAM files
# -@ number of threads to use (not necessary)
# -b output bam file
# -q 20: keep reads with mapping quality above 20
# -f 0x0002 only keep proper pairs (as defined by bwa)
# -F 0x0004 remove reads that are not mapped
# -F 0x0008 remove reads with an un-mapped mate
# -F 0x100 remove secondary alignments
samtools view -b -q 20 \
    -f 0x0002 \
    -F 0x0004 \
    -F 0x0008 \
    -F 0x100 \
    ${OUTPUT}.markdup.bam > ${OUTPUT}.markdup.Q20.bam

# gather statistics
samtools flagstat ${OUTPUT}.markdup.Q20.bam > ${OUTPUT}.markdup.Q20.flagstat