#!/bin/bash

REF=~/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa.gz

cd ~
bcftools mpileup -a AD,DP,SP -Ou -f $REF \
./results/align/Plomion/*_sort.bam | bcftools call -f GQ,GP \
-mO z -o ./results/VCF/quercus_bcftools.vcf.gz
