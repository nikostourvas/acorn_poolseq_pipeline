#!/bin/bash

#This script takes chunks of a whole-genome vcf that are placed together in a directory. They are merged and transformed to a set of biallelic SNPs from which all INDEL-adjacent sites are removed. 

UNMERGED_VCF_DIR=$1 #The path to the directory in which all the unmerged VCF files are stored
VCF_PREFIX=$2 #The name for the merged VCF. Prefix only, no need for ".vcf". 
THREADS=$3 #The number of CPU cores to use

# Merge the multiple small VCF files that were produced in the previous step
# into one large VCF
bash /mnt/acorn_poolseq_pipeline/Concatenate_VCFs.sh ${UNMERGED_VCF_DIR} ${VCF_PREFIX} ${THREADS} 

# Filter SNPs detected close to INDELs, as suggested by VarScan manual
bash /mnt/acorn_poolseq_pipeline/snp_indel_rm.sh ${UNMERGED_VCF_DIR}/${VCF_PREFIX}_SNP.vcf

# Filter out multi-allelic SNPs
bash /mnt/acorn_poolseq_pipeline/vcf_biallelic.sh ${UNMERGED_VCF_DIR}/${VCF_PREFIX}IndelFilteredSNPs.vcf
