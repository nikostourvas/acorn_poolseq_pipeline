#!/bin/bash

#Lars Littmann
#21 Nov 2023
#Randomly subsample a specified number of SNPs from a VCF file. The output file is placed in the input file's directory.
#This script needs to be run in the GATK container

#INPUTS
INPUT_VCF=$1 #UNZIPPED (!) VCF from which to take a random sample of SNPs
SUBSAMPLE_RATIO=$2

#PREPARE OUTPUT
OUTPUT=${INPUT_VCF/.vcf/}_Subsampled_${SUBSAMPLE_RATIO}.vcf

gatk SelectVariants -R ./REFERENCE/Qrob_PM1N.fa -V ${INPUT_VCF} -O ${OUTPUT} --select-random-fraction ${SUBSAMPLE_RATIO}
