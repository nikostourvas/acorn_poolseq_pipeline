#!/bin/bash

# declare variables
RESULTS=../results/VCF
REF=../reference/Qrob_PM1N.fa

# Remove SNPs close to InDels & perform further SNP filtering
java -jar /usr/share/java/varscan.jar filter $RESULTS/Qrob_total.snp.vcf \
    --min-var-freq 0.025 --p-value 0.05 --min-avg-qual 20 \
    --min-coverage 30 --min-reads2 1 \
    --indel-file $RESULTS/Qrob_total.indel.vcf \
    --output-file $RESULTS/Qrob_total_filter.snp.vcf \
    2> $RESULTS/varscan_SNP_filter.err

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