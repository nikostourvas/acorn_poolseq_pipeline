#!/bin/bash

# fastqc
bash fastqc.sh untrimmed_fastq data

# Trimmomatic
parallel --verbose -j 2 \
	'bash trimmomatic_parallel.sh {}' :::: inds

# fastqc again
bash fastqc.sh trimmed_fastq results

# Map
bwa index /home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa
samtools faidx /home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa

parallel --verbose -j 2 \
	'bash mapping.sh {}' :::: inds

parallel --verbose -j 2 \
	'bash bam_cleaning.sh {}' :::: inds

# Variant calling
bash variant_calling_samtools-VarScan.sh
