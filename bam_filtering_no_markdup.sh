#!/bin/bash

# Nikos Tourvas & Lars Littmann
# 2023
# Filter raw BAM files with samtools view
# This version of the script does NOT remove duplicates!

# declare variables
IND=$1 #The individual that is processed
BAM=/data/genetics_tmp/results/mapped_reads/${IND} #The location of the bam files, by individual.

# Quick explanation of what individual commands do:
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
samtools sort -u ${BAM}.raw.bam 2> ${BAM}.sort.err \
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

##################################################################
#A handy way to run this script is as follows:
#Create a .txt file that contains all the sample names of samples you wish to process. These sample names have to match with your file naming. 
#Run the following parallel command. It initiates this bash script for each sample. As soon as one script finishes running, it starts the next until it has cycled through all file names in the .txt file.
#These commands are not part of the script, but can be copied into the command line to run the script. 

#parallel --verbose -j 90 \
#	'bash bam_filtering_no_markdup.sh {}' :::: /data/genetics_tmp/results/fastp_dedup_trim/samplenames_bySize.txt
