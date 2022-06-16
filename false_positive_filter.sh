#!/bin/bash

# declare variables
BAM_FILE="/home/tourvasn/ngs_training/results/align/54977_ID1593_8-209-Tube_S1_L001.sort.Q20.markdup.bam"
RESULTS="/home/tourvasn/ngs_training/results/VCF"
REF="/home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa"
BASE="/home/tourvasn/ngs_training/results/VCF/"
SNPVCFS=$RESULTS/*snp.vcf
INDELVCFS=$RESULTS/*indel.vcf
FPFILTER="/home/tourvasn/other_software"

# Obtain metrics for the list of variants
#bam-readcount -w 2 -q 1 -b 20 -f $REF \
#    $BAM_FILE \
#    >$RESULTS/varScan.variants.readcounts

# Run the FPfilter accessory script
perl $FPFILTER/fpfilter.pl $RESULTS/Qrob_total_filter.snp.vcf \
    $RESULTS/varScan.variants.readcounts \
    --output-basename Qrob_total_fp_filter