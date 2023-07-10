#!/bin/bash

# declare variables
RESULTS=../results/VCF
REF=/mnt/data/reference/Qrob_PM1N.fa

# keep only biallelic SNPs
bcftools view --threads 32 -m2 -M2 -v snps $RESULTS/Qrob_total_filter.snp.vcf \
    -o $RESULTS/Qrob_total_filter2.snp.vcf

# produce statistics
bcftools stats --threads 32 --fasta-ref $REF $RESULTS/Qrob_total_filter2.snp.vcf \
    > $RESULTS/bcftools_stats.txt