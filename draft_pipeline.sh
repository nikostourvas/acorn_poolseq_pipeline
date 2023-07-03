#!/bin/bash

# Check FASTQ file integrity
# Generate a hash for each file
md5sum ../data/untrimmed_fastq/*.fastq.gz >  ../data/untrimmed_fastq/hashList.txt &&

# fastqc - Quality control for raw sequences
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh untrimmed_fastq data &&

# make list of samples to allow for parallel analysis of multiple samples
for i in ../data/untrimmed_fastq/*_R1_001.fastq.gz
do
	echo $(basename ${i%_R*})
done > ../data/inds &&

# Trimmomatic - filter out bad reads or bad parts of a read
# You can change number of cores per job by accessing trimmomatic_parallel.sh script.
# For example you can set trimmomatic to use 4 cores. 
# ATTN: So if you ask parallel to run 2 jobs concurrently, you need 2 x 4 = 8 cores.
parallel --verbose -j 2 \
	'bash trimmomatic_parallel.sh {}' :::: ../data/inds &&

# fastqc again - Quality control for trimmed sequences
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh trimmed_fastq results &&

# Create a unified report with MultiQC
multiqc ../results -o ../results &&

# Map reads to reference genome
# You can change number of cores per job by accessing mapping.sh script.
bwa index ../data/reference/Qrob_PM1N.fa &&
samtools faidx ../data/reference/Qrob_PM1N.fa &&

parallel --verbose -j 2 \
	'bash mapping.sh {}' :::: ../data/inds &&

parallel --verbose -j 2 \
	'bash bam_cleaning.sh {}' :::: ../data/inds &&

# make list of all chromosomes & scaffolds for next step
grep 'Qrob' ../data/reference/Qrob_PM1N.fa | \
 cut -c2- > ../data/regions && 

# Variant calling
parallel --verbose -j 20 \
	'bash variant_calling_parallel.sh {}' :::: ../data/regions &&

# Merge the multiple small VCF files that were produced in the previous step
# into one large VCF
bash vcf_merge.sh &&

# Filter SNPs detected close to INDELs, as suggested by VarScan manual
bash snp_indel_rm.sh &&

# Filter out multi-allelic SNPs
bash vcf_biallelic.sh &&

# Create a unified final report with MultiQC
multiqc --force ../results -o ../results