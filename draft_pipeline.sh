#!/bin/bash

# Check FASTQ file integrity
# Generate a hash for each file
mkdir -p ../results
md5sum /mnt/data/untrimmed_fastq/*.fastq.gz >  ../results/hashList.txt &&

# fastqc - Quality control for raw sequences
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh untrimmed_fastq /mnt/data &&

# make list of samples to allow for parallel analysis of multiple samples
for i in /mnt/data/untrimmed_fastq/*_R1_001.fastq.gz
do
	echo $(basename ${i%_R*})
done > /mnt/data/inds &&

# Trimmomatic - filter out bad reads or bad parts of a read
# You can change number of cores per job by accessing trimmomatic_parallel.sh script.
# For example you can set trimmomatic to use 4 cores. 
# ATTN: So if you ask parallel to run 2 jobs concurrently, you need 2 x 4 = 8 cores.
parallel --verbose -j 2 \
	'bash trimmomatic_parallel.sh {} 1>/mnt/results/trimmed_fastq_Batch1/logs/{}.out 2> /mnt/results/trimmed_fastq_Batch1/logs/{}.error' :::: /mnt/AcornSeqdata/inds_Batch1.txt &&

# fastqc again - Quality control for trimmed sequences
# You can change number of cores per job by accessing fastqc.sh script.
bash fastqc.sh trimmed_fastq ../results &&

# Create a unified report with MultiQC
multiqc ../results -o ../results &&

# Map reads to reference genome
# You can change number of cores per job by accessing mapping.sh script.
bwa index /mnt/data/reference/Qrob_PM1N.fa &&
samtools faidx /mnt/data/reference/Qrob_PM1N.fa &&

parallel --verbose -j 12 \
	'bash map.sh {}' :::: ../AcornSeqdata/inds_Batch1.txt

parallel --verbose -j 7 \
	'bash bam_filtering.sh {}' :::: ../AcornSeqdata/inds_Batch1.txt

# make list of all chromosomes & scaffolds for next step
grep 'Qrob' /mnt/data/reference/Qrob_PM1N.fa | \
 cut -c2- > /mnt/data/regions && 

# Variant calling
parallel --verbose -j 4 \
	'bash variant_calling_parallel.sh {}' :::: /mnt/data/regions &&

# Merge the multiple small VCF files that were produced in the previous step
# into one large VCF
bash vcf_merge.sh &&

# Filter SNPs detected close to INDELs, as suggested by VarScan manual
bash snp_indel_rm.sh &&

# Filter out multi-allelic SNPs
bash vcf_biallelic.sh &&

# Create a unified final report with MultiQC
multiqc --force ../results -o ../results
