#!/bin/bash

VCF=../results/VCF/Qrob_total_filter2.snp.vcf
OUTPUT=../results/VCF/acorn_poolseq
MISSING=0.1

# Use plink1.9 (i)to convert VCF to bed and
#              (ii)to remove missing data (--geno argument) 
# The argument --double-id --allow-extra-chr is necessary because we are
# working non-human sequences
# The argument --bp-space can also be used to exclude one variant from 
# each pair closer than the give bp count (similarly to vcftools --thin)
# e.g. Leroy et al. 2020 kept one allele every 15kb for ABC analyses
plink1.9 --vcf $VCF \
         --double-id --allow-extra-chr \
         --set-missing-var-ids @_# \
         --geno $MISSING \
         --make-bed \
         --out $OUTPUT