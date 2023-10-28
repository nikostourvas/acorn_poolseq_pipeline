#!/bin/bash

#Concatenate VCFs
#Lars Littmann
#2023-09-04
#Takes (gz)vcfs that contain the same samples and combines them into one file. In the ACORN project, this applies for VCFs that have been generated for
#particular chunks of the reference genome and now require merging. Operates both on the indel and on the snp vcf files.

UNMERGED_VCF_DIR=$1 #The path to the directory in which all the unmerged VCF files are stored
VCF_PREFIX=$2 #The name for the merged VCF
THREADS=$3 #The number of CPU cores to use

#Create a list of filepaths containing each of the VCFs we wish to merge. One for SNPs, one for Indels
realpath ${UNMERGED_VCF_DIR}/chunk*.varScan.snp.vcf.gz | sort -V > ${UNMERGED_VCF_DIR}/SNP_VCF_List.txt
realpath ${UNMERGED_VCF_DIR}/chunk*.varScan.indel.vcf.gz | sort -V > ${UNMERGED_VCF_DIR}/INDEL_VCF_List.txt

#Concatenate the VCFs containing SNPs using bcftools.
#This implementation does not allow for overlaps
#-f takes input from the list of filepaths we created earlier.
#-a allows for discontiguous regions to be accepted (manual check: This only applies to the very edges of chromosomes. No harm there).
#-O z sets the output type to gzipped vcf, and v to uncompressed vcf
#-o is the output prefix (specified by user)
bcftools concat -f ${UNMERGED_VCF_DIR}/SNP_VCF_List.txt \
    -a -O v -o ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_SNP.vcf \
    --threads ${THREADS} 2> ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_SNP.vcf.err

#Concatenate the VCFs containing INDELs using bcftools.
bcftools concat -f ${UNMERGED_VCF_DIR}/INDEL_VCF_List.txt \
    -a -O v -o ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_INDEL.vcf \
    --threads ${THREADS} 2> ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_INDEL.vcf.err

#Remove some intermediate files
rm ${UNMERGED_VCF_DIR}/SNP_VCF_List.txt
rm ${UNMERGED_VCF_DIR}/INDEL_VCF_List.txt
