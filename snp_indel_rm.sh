#!/bin/bash

#Nikos Tourvas
#2023
#Script removes SNPs from a VCF that are in close proximity to INDEL regions. The script takes in two VCFs; one that contains SNPs and one that contains INDELs.
#Some additional SNP filtering is also performed in this script.
#This script can be run by itself, or as part of the script post_variant_calling.sh
#https://github.com/nikostourvas/acorn_poolseq_pipeline/blob/singularity/post_variant_calling.sh

# declare variables
SNP_VCF=$1 #Path to a vcf containing SNPs (tested on VarScan VCFs)
INDEL_VCF=${SNP_VCF/_SNP.vcf/_INDEL.vcf} #Automatically find a VCF containing INDELs. This file should be in the same directory as the SNP vcf. 
OUTDIR=$(dirname ${SNP_VCF}) #Specify an output directory. 
REF=/mnt/reference/Qrob_PM1N.fa #The path to a reference sequence. In this case the Plomion et al (2018) Quercus robur reference genome with mitochondrial and chloroplast genomes added. 

# Remove SNPs close to InDels & perform further SNP filtering
# Remove SNPs that have been called in sites with low read quality (avg. Phred score <20).
# Implement a minimum read depth threshold of 20 (in our case, might not be applicable for other projects).
# Remove SNPs called with a p-value above 0.05. 
# Remove SNPs for which the alt has fewer than one reads (practically this filter is redundant. Merely to make this fact explicit).
java -Xmx128g -jar /usr/share/java/varscan.jar filter ${SNP_VCF} \
    --min-var-freq 0.00 --p-value 0.05 --min-avg-qual 20 \
    --min-coverage 20 --min-reads2 1 \
    --indel-file ${INDEL_VCF} \
    --output-file ${OUTDIR}/$(basename ${SNP_VCF/_SNP.vcf/IndelFilteredSNPs.vcf}) \
    2> ${OUTDIR}/varscan_SNP_filter.err

# OPTIONAL TODO!
# Use VarScan false positive filter
# The scientific basis of this filter is described in the VarScan 2 publication. It will improve
#the precision of variant and mutation calling by removing artifacts associated with short-read alignment.
#-For somatic mutations, generate bam-readcounts with the Tumor BAM. For LOH and Germline, generate readcounts with the Normal BAM
#-For de novo mutations (trio calling), generate readcounts with the child BAM.
# The filter requires the bam-readcount utility: https://github.com/genome/bam-readcount
# Individual BAM files need to be merged into a single large BAM file!!!

# Merge bams
#samtools merge -r $RESULTS/all.bam $BAM_FILES \
        #--threads 4
# Readcount BAM
#bam-readcount -w 1 -q 1 -b 20 -f $REF $RESULTS/all.bam > $RESULTS/readcount.tsv

#java -jar /usr/share/java/varscan.jar fpfilter $RESULTS/Qrob_total_filter.snp.vcf \
#     $RESULTS/readcount.tsv --output-file $RESULTS/Qrob_total_filterfp.snp.vcf
