#!/bin/bash

# declare variables
IND=$1 #The individual that is processed
BAM=/data/genetics_tmp/results/mapped_reads/${IND} #The location of the bam files, by individual.

# Filter raw BAM files with samtools view
# This version of the script does NOT remove duplicates!

# Quick explanation of what individual commands do:
# collate: groups the reads by read name in the bam file. This puts the read 
# pairs close together so that fixmate can work 
# fixmate: fix read-mates so that they both have the same sets of attributes for
# the subsequent preprocessing
# sort: sort reads again based on genomic coordinates and
# Then produce filtered BAM files with samtools view and the following arguments
# -@ number of threads to use (not necessary)
# -b output bam file
# -L keep only contigs specified in the provided bed file (contigs >500 kb)
# -q 20: keep reads with mapping quality above 20
# -f 0x0002 only keep proper pairs (as defined by bwa)
# -F 0x0004 remove reads that are not mapped
# -F 0x0008 remove reads with an un-mapped mate
# -F 0x100 remove secondary alignments
samtools collate -O -u ${BAM}.raw.bam 2> ${BAM}.collate.err \
    | samtools fixmate -m -u - - 2> ${BAM}.fixmate.err \
    | samtools sort -u - 2> ${BAM}.sort.err \
    | samtools view -h -b \
        -L /data/genetics_tmp/REFERENCE/contigsover500kb.bed \
        -q 20 \
        -f 0x0002 \
        -F 0x0004 \
        -F 0x0008 \
        -F 0x100 \
        -o ${BAM}.filtered.bam - 2> ${BAM}.filter-view.err

# index bam files
samtools index -@ 1 ${BAM}.filtered.bam

# gather statistics
# -@ number of cores
samtools flagstat -@ 1 ${BAM}.filtered.bam > ${BAM}.filtered.flagstat