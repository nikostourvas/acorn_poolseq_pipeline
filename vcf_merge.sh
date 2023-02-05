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
# TODO: multithreading on bcftools index
for SNPVCF in ${SNPVCFS[@]}
do
    bcftools index ${SNPVCF}.gz
done

for INDELVCF in ${INDELVCFS[@]}
do
    bcftools index ${INDELVCF}.gz
done

# Merge VCFs
# TODO: multithreading on bcftools concat
bcftools concat $RESULTS/Qrob*.snp.vcf.gz -o $RESULTS/Qrob_total.snp.vcf
bcftools concat $RESULTS/Qrob*.indel.vcf.gz -o $RESULTS/Qrob_total.indel.vcf

# Filter SNPs close to InDels
java -jar /usr/share/java/varscan.jar filter $RESULTS/Qrob_total.snp.vcf \
    --min-var-freq 0.0055 --p-value 0.05 --min-avg-qual 20 \
    --min-coverage 50 --min-reads2 2 \
    --indel-file $RESULTS/Qrob_total.indel.vcf \
    --output-file $RESULTS/Qrob_total_filter.snp.vcf

# keep only biallelic SNPs
# TODO: multithreading on bcftools view possible???
bcftools view -m2 -M2 -v snps $RESULTS/Qrob_total_filter.snp.vcf \
    -o $RESULTS/Qrob_total_filter2.snp.vcf