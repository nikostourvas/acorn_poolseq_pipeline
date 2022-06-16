#!/bin/bash

# declare variables
REF=/home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa
BAM=/home/tourvasn/ngs_training/results/align/54977_ID1593_8-209-Tube_S1_L001
RES=/home/tourvasn/ngs_training/results/VCF/209.vcf

java -Xmx8g -jar /programs/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
    HaplotypeCaller \
    -R ${REF} \
    -I ${BAM}.sort.Q20.markdup.bam \
    --sample-ploidy 88 \
    --max-genotype-count 3 \
    -O ${RES}