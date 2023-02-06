#!/bin/bash

# declare variables
RESULTS=../results/VCF

# keep only biallelic SNPs
bcftools view --threads 32 -m2 -M2 -v snps $RESULTS/Qrob_total_filter.snp.vcf \
    -o $RESULTS/Qrob_total_filter2.snp.vcf