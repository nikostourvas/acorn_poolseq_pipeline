#!/bin/bash

# declare variables
RESULTS=../results/VCF

# Filter SNPs close to InDels
java -jar /usr/share/java/varscan.jar filter $RESULTS/Qrob_total.snp.vcf \
    --min-var-freq 0.0055 --p-value 0.05 --min-avg-qual 20 \
    --min-coverage 50 --min-reads2 2 \
    --indel-file $RESULTS/Qrob_total.indel.vcf \
    --output-file $RESULTS/Qrob_total_filter.snp.vcf