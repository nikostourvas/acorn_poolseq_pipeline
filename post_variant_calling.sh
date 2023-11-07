#!/bin/bash

#UNMERGED_VCF_DIR=$1 #The path to the directory in which all the unmerged VCF files are stored
#VCF_PREFIX=$2 #The name for the merged VCF
#THREADS=$3 #The number of CPU cores to use

# Merge the multiple small VCF files that were produced in the previous step
# into one large VCF
bash /mnt/acorn_poolseq_pipeline/Concatenate_VCFs.sh ${1} ${2} ${3} 

# Filter SNPs detected close to INDELs, as suggested by VarScan manual
bash /mnt/acorn_poolseq_pipeline/snp_indel_rm.sh ${1}/${2}_SNP.vcf

# Filter out multi-allelic SNPs
bash /mnt/acorn_poolseq_pipeline/vcf_biallelic.sh ${1}/${2}IndelFilteredSNPs.vcf