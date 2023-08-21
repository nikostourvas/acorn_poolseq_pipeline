#!/bin/bash

# declare variables
RESULTS=../results/VCF
REF=../reference/Qrob_PM1N.fa
THREADS=12

# keep only biallelic SNPs
bcftools view --threads ${THREADS} -m2 -M2 -v snps \
    $RESULTS/Qrob_total_filter.snp.vcf \
    -Oz -o $RESULTS/Qrob_total_filter2.snp.vcf.gz \
    2> ${RESULTS}.bcftools.biallelic.vcf.err

# produce statistics
bcftools stats --threads ${THREADS} --fasta-ref $REF \
    $RESULTS/Qrob_total_filter2.snp.vcf.gz \
    > $RESULTS/bcftools_stats.txt