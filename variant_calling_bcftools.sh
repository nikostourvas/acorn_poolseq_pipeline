#!/bin/bash

REF=~/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa.gz

cd ~
bcftools mpileup -a AD,DP,SP -Ou -f $REF \
./ngs_training/results/align/Plomion/*_sort.bam | bcftools call -f GQ,GP \
-mO z -o ./ngs_training/results/VCF/quercus_bcftools.vcf.gz
