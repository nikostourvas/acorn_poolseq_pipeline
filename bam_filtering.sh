#!/bin/bash

# declare variables
IND=$1 #The individual that is processed
BAM=../results/align_Batch1/${IND} #The location of the bam files, by individual.

# Remove PCR and optical duplicates with samtools markdup
# This is the recommended way to mimic Picard tools MarkDuplicates with samtools
# See https://www.htslib.org/algorithms/duplicate.html
# After that the BAM file is filtered with samtools view

# Quick explanation of what individual commands do:
# collate: groups the reads by read name in the bam file. This puts the read 
# pairs close together so that fixmate can work 
# fixmate: fix read-mates so that they both have the same sets of attributes for
# the subsequent preprocessing
# sort: sort reads again based on genomic coordinates and
# markdup: remove duplicates (-r argument)
# Then produce filtered BAM files with samtools view and the following arguments
# -@ number of threads to use (not necessary)
# -b output bam file
# -q 20: keep reads with mapping quality above 20
# -f 0x0002 only keep proper pairs (as defined by bwa)
# -F 0x0004 remove reads that are not mapped
# -F 0x0008 remove reads with an un-mapped mate
# -F 0x100 remove secondary alignments
samtools collate -O -u ${BAM}.raw.bam 2> ${BAM}.collate.err \
    | samtools fixmate -m -u - - 2> ${BAM}.fixmate.err \
    | samtools sort -u - 2> ${BAM}.sort.err\
    | samtools markdup -u -r -f ${BAM}_markdup_stats_file.txt \
        -S -d 2500 --mode s --include-fails - - 2> ${BAM}.markdup.err \
    | samtools view -b -q 20 \
        -f 0x0002 \
        -F 0x0004 \
        -F 0x0008 \
        -F 0x100 \
        -o ${BAM}.markdup.Q20.bam - 2> ${BAM}.markdup-view.err

# gather statistics
# -@ number of cores
samtools flagstat -@ 1 ${BAM}.markdup.Q20.bam > ${BAM}.markdup.Q20.flagstat

# index bam files
samtools index -@ 1 ${BAM}.markdup.Q20.bam
