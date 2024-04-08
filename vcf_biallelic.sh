#!/bin/bash

#Nikos Tourvas & Lars Littmann
#2023
#Script takes any vcf and retains only bi-allelic SNPs in an output VCF. The script also outputs statistics for the newly created vcf.
#This script can be run by itself, or as part of the script post_variant_calling.sh
#https://github.com/nikostourvas/acorn_poolseq_pipeline/blob/singularity/post_variant_calling.sh

# declare variables
INPUT_VCF=$1
REF=/mnt/reference/Qrob_PM1N.fa #The path to a reference sequence. In this case the Plomion et al (2018) Quercus robur reference genome with mitochondrial and chloroplast genomes added. 
THREADS=60
OUTDIR=$(dirname ${INPUT_VCF})

# keep only biallelic SNPs using bcftools and setting both the minimum (-m) and maximum (-M) number of variants at a site at 2.
bcftools view --threads ${THREADS} -m2 -M2 -v snps \
    ${INPUT_VCF} \
    -Oz -o ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_Biallelic.vcf.gz}) \
    2> ${OUTDIR}/bcftools_biallelic_vcf.err

# produce statistics for the newly created biallelic VCF. 
bcftools stats --threads ${THREADS} --fasta-ref ${REF} \
    ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_Biallelic.vcf.gz}) \
    > ${OUTDIR}/$(basename ${INPUT_VCF/.vcf/_bcftools_stats.txt})
