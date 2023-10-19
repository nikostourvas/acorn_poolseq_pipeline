#!/bin/bash

# declare variables
INPUT_VCF=$1
REF=../reference/Qrob_PM1N.fa
THREADS=12
OUTDIR=$(dirname ${INPUT_VCF})

# keep only biallelic SNPs
bcftools view --threads ${THREADS} -m2 -M2 -v snps \
    ${INPUT_VCF} \
    -Oz -o ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_Biallelic.vcf.gz}) \
    2> ${OUTDIR}/bcftools_biallelic_vcf.err

# produce statistics
bcftools stats --threads ${THREADS} --fasta-ref ${REF} \
    ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_Biallelic.vcf.gz}) \
    > ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_bcftools_stats.txt})
