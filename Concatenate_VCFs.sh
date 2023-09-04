#!/bin/bash

#Concatenate VCFs
#Lars Littmann
#2023-09-04
#Takes (gz)vcfs that contain the same samples and combines them into one file. In the ACORN project, this applies for VCFs that have been generated for
#particular chunks of the reference genome and now require merging. Operates both on the indel and on the snp vcf files.

UNMERGED_VCF_DIR=$1 #The path to the directory in which all the unmerged VCF files are stored
VCF_PREFIX=$2 #The name for the merged VCF

#Create a list of filepaths containing each of the VCFs we wish to merge. One for SNPs, one for Indels
realpath ${UNMERGED_VCF_DIR}/chunk*.varScan.snp.vcf.gz > ${UNMERGED_VCF_DIR}/SNP_VCF_List.txt
realpath ${UNMERGED_VCF_DIR}/chunk*.varScan.indel.vcf.gz > ${UNMERGED_VCF_DIR}/INDEL_VCF_List.txt

#Concatenate the VCFs containing SNPs using bcftools
bcftools concat -f ${UNMERGED_VCF_DIR}/SNP_VCF_List.txt -O z -o ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_SNP.vcf.gz

#Concatenate the VCFs containing INDELs using bcftools
bcftools concat -f ${UNMERGED_VCF_DIR}/INDEL_VCF_List.txt -O z -o ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_INDEL.vcf.gz
