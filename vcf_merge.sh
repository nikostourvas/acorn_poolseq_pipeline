#!/bin/bash

# declare variables
RESULTS=../results/VCF
SNPVCFS=$RESULTS/*snp.vcf
INDELVCFS=$RESULTS/*indel.vcf
THREADS=12

# Merge VCF files containing SNPs from each genomic region
bcftools concat --threads ${THREADS} $RESULTS/Qrob*.snp.vcf.gz \
    -o $RESULTS/Qrob_total.snp.vcf \
    2> ${RESULTS}/${REGION}.bcftools.snp.merge.vcf.err

# Merge VCF files containing INDELs from each genomic region
bcftools concat --threads ${THREADS} $RESULTS/Qrob*.indel.vcf.gz \
    -o $RESULTS/Qrob_total.indel.vcf \
    2> ${RESULTS}/${REGION}.bcftools.indel.merge.vcf.err

# Here it would be optimal to output compressed VCF files. However the
# downstream 'varscan filter' command does not work with compressed files.