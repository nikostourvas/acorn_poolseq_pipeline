#!/bin/bash

REF=~/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa

cd ~
bcftools mpileup -a AD,DP,SP -Ou -f $REF --threads 20\
./ngs_training/results/align/Plomion/*_sort.bam | bcftools call -f GQ,GP --ploidy 40 \
--threads 20 \
-mO z -o ./ngs_training/results/VCF/quercus_bcftools.vcf.gz
