#!/bin/bash

#Lars Littmann 
#12/02/2024
#Basic filtering to be applied to all PoolSeq VCFs created using VarScan for the ACORN project.
#This script does not remove SNPs situated near InDels. This is done by an upstream script. 
#This script takes the following inputs: 
#A gzipped PoolSeq VCF. 
#Threshold value for minimum read depth
#Threshold value for maximum allowed missingness across samples (pools)

#The script calculates and outputs a maximum read depth score. This is the mean plus 2 times the standard deviation from the mean. 

#Read inputs from command line

VCF_IN=$1
MIN_RD=$2
MISSINGNESS=$3

#Specify further variables

FILENAME=$(basename ${VCF_IN/.vcf.gz/})
OUT_DIR=$(dirname ${VCF_IN})
INT_DIR=${OUT_DIR}/IntermediateFiles_${FILENAME}
MISSINGNESS_STRING=$(echo "'F_MISSING>${MISSINGNESS}'") #This is necessary because bcftools does not allow you to place a variable between quotes, but echo does.

#Create intermediate directory

mkdir -p ${INT_DIR}

#Change the version 4.3 vcf produced by VarScan to a 4.2 vcf.
zcat ${VCF_IN} | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/g' > ${INT_DIR}/${FILENAME}_V42.vcf

#Extract the mean depth across all pools per site using vcftools
vcftools --vcf ${INT_DIR}/${FILENAME}_V42.vcf --site-mean-depth --out ${INT_DIR}/ReadDepths

#use awk to get the mean and standard deviation in mean site read depth.
MEAN=$(tail -n +2 ${INT_DIR}/ReadDepths.ldepth.mean | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
SD=$(tail -n +2 ${INT_DIR}/ReadDepths.ldepth.mean | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)*2)}')

#Use the above values to calculate a threshold for maximum read depth filtering
MAX_RD_FLOAT=$((${MEAN}+2*${SD}))
#This is not quite enough. bcftools can only take integers, so we need to round this number.
MAX_RD_INT= $(echo ${MAX_RD_FLOAT} | awk '{print int($1+0.5)}')

#Report the values in an intermediate file
echo "The average read depth calculated using vcftools --site-mean-depth and the standard deviation are as follows:" > ReadDepthSummary.txt
echo "Mean = ${MEAN}, SD = ${SD}" >> ReadDepthSummary.txt
echo "Making the maximum read depth threshold (Mean+2xSD) = ${MAX_RD_FLOAT}, which is rounded to ${MAX_RD_INT}" >> ReadDepthSummary.txt

#Now that we have all the values we want, let's start filtering. First up is minimum and maximum read depth. We filter for this using VCFtools
vcftools --vcf ${INT_DIR}/${FILENAME}_V42.vcf --minDP ${MIN_RD} --maxDP ${MAX_RD_INT} \
--recode --recode-INFO-all --out ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxDP${MAX_RD_INT}

#Next up, we filter out sites with poor quality (which won't filter out much), and filter for missingness using bcftools
bcftools view -i 'MIN(FMT/GQ)>15' ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxDP${MAX_RD_INT}.recode.vcf | \
bcftools filter -e ${MISSINGNESS_STRING} > ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxDP${MAX_RD_INT}_Missingness${MISSINGNESS}.vcf


#As I was writing this script, I got insanely distracted at some point and made this chicken. Enjoy.

#   _\/_
#  / *  \      ////
# <      \____/  /   buh-gok
#  \            /     
#   \__     ___/     _
#      \___/        / \
#        |         |   |
#       //          \_/
