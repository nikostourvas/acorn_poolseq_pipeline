#!/bin/bash

# fastqc
bash fastqc.sh untrimmed_fastq data

# make list of samples
for i in ../data/untrimmed_fastq/*_R1_001.fastq.gz
do
	echo $(basename ${i%_R*})
done > ../data/inds

# Trimmomatic
parallel --verbose -j 2 \
	'bash trimmomatic_parallel.sh {}' :::: ../data/inds

# fastqc again
bash fastqc.sh trimmed_fastq results

# Map
bwa index ../data/reference/Qrob_PM1N.fa
samtools faidx ../data/reference/Qrob_PM1N.fa

parallel --verbose -j 2 \
	'bash mapping.sh {}' :::: ../data/inds

parallel --verbose -j 2 \
	'bash bam_cleaning.sh {}' :::: ../data/inds

# make list of all chromosomes & scaffolds
grep 'Qrob' ../data/reference/Qrob_PM1N.fa | \
 cut -c2- > ../data/regions 

# Variant calling
#####REMOVE IF NOT USEDbash variant_calling_samtools-VarScan.sh
parallel --verbose -j 20 \
	'bash variant_calling_parallel.sh {}' :::: ../data/regions

bash vcf_merge.sh

# filter for paralogs
# extremely RAM hungry - DO NOT RUN
#singularity exec ~/singularity_images/bam-readcount.sif false_positive_filter.sh