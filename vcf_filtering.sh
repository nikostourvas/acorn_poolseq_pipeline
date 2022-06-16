#!/bin/bash

# declare variables
BASE=/home/tourvasn/ngs_training/results_backup/VCF
REF=/home/tourvasn/ngs_training/data/reference/Qrob_PM1N.fa
INPUT=$BASE/Qrob_209-220.varScan.snp.vcf.gz
SUBSET_VCF=$BASE/Qrob_subset.vcf
OUT=$BASE/Qrob_subset.table

# zip the vcf
#bgzip -k $BASE/Qrob_209-220.varScan.snp.vcf
#index vcf
#bcftools index $INPUT
#tabix $INPUT

# Randomly subsampling VCF
#bcftools view $INPUT | \
#    vcfrandomsample -r 0.025 > $SUBSET_VCF
#    bgzip > $SUBSET_VCF

# calculate allele frequency
#vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2
# calculate mean depth per pool
#vcftools --gzvcf $SUBSET_VCF --depth --out $OUT
# calculate mean depth per site
#vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT
# calculate proportion of missing data per pool - does this make sense for pool?
#vcftools --gzvcf $SUBSET_VCF --missing-pool --out $OUT
# calculate proportion of missing data per site
#vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# create dict reference genome file for gatk
#samtools dict $REF > /home/tourvasn/ngs_training/data/reference/Qrob_PM1N.dict

java -Xmx4g -jar /programs/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar \
    VariantsToTable \
    -V $SUBSET_VCF \
    -F CHROM -F POS -F TYPE -F REF -F ALT -GF GT -GF GQ -GF FREQ \
    -O $OUT --split-multi-allelic