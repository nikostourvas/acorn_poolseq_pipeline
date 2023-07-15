#!/bin/bash

# declare variables
RESULTS=../results/VCF
SNPVCFS=$RESULTS/*snp.vcf
INDELVCFS=$RESULTS/*indel.vcf

# some scaffolds will have no snps/indels
# delete vcf files originating from these scaffolds
find $RESULTS -maxdepth 1 -type f -empty -print -delete

# zip VCFs
for SNPVCF in ${SNPVCFS[@]}
do
    bgzip -k $SNPVCF
done

for INDELVCF in ${INDELVCFS[@]}
do
    bgzip -k $INDELVCF
done

# Index VCFs
for SNPVCF in ${SNPVCFS[@]}
do
    bcftools index --threads 4 ${SNPVCF}.gz
done

for INDELVCF in ${INDELVCFS[@]}
do
    bcftools index --threads 4 ${INDELVCF}.gz
done

# Merge VCFs
bcftools concat --threads 4 $RESULTS/Qrob*.snp.vcf.gz -o $RESULTS/Qrob_total.snp.vcf
bcftools concat --threads 4 $RESULTS/Qrob*.indel.vcf.gz -o $RESULTS/Qrob_total.indel.vcf