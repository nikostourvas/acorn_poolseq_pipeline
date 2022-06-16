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

# make list of all chromosomes & scaffolds
grep 'Qrob' /home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa | \
 cut -c2- > regions 

# Variant calling
bash variant_calling_samtools-VarScan.sh
parallel --verbose -j 20 \
	'bash variant_calling_parallel.sh {}' :::: regions

bash vcf_merge.sh

# filter for paralogs
# extremely RAM hungry - DO NOT RUN
#singularity exec ~/singularity_images/bam-readcount.sif false_positive_filter.sh