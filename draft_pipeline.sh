#!/bin/bash

# fastqc
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh untrimmed_fastq data

# make list of samples
for i in ../data/untrimmed_fastq/*_R1_001.fastq.gz
do
	echo $(basename ${i%_R*})
done > ../data/inds

# Trimmomatic
# You can change number of cores per job by accessing trimmomatic_parallel.sh script.
# For example you can set trimmomatic to use 4 cores. 
# ATTN: So if you ask parallel to run 2 jobs concurrently, you need 2 x 4 = 8 cores.
parallel --verbose -j 2 \
	'bash trimmomatic_parallel.sh {}' :::: ../data/inds

# fastqc again
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh trimmed_fastq results

# Map
# You can change number of cores per job by accessing mapping.sh script.
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

bash snp_indel_rm.sh

bash vcf_biallelic.sh

# filter for paralogs
# extremely RAM hungry - DO NOT RUN
#singularity exec ~/singularity_images/bam-readcount.sif false_positive_filter.sh