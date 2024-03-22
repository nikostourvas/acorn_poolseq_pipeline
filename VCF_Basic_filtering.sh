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

VCF_IN=$1 #VCF can be either gzipped or full-fat
MIN_RD=$2 #The minimum read depth. Integer
MISSINGNESS=$3 # The maximum missingness allowed for a site, expressed as a decimal (example 0.1)
MAF=$4 #The minimum allele frequency for a SNP to be included, expressed as a decimal (example 0.05)

#Specify further variables

HD_MASK=/data/genetics_tmp/REFERENCE/MASKS/HDplot_Mask_D10_H06_NoWindow.txt #Specify the location of a relevant mask generated using HDplot. Requires manual editing to change the mask.
FILENAME=$(basename ${VCF_IN/.vcf.gz/}) #The base for output file names.
OUT_DIR=$(dirname ${VCF_IN}) #Specify where to output. In this case, it is the same directory in which our vcf is.
MISSINGNESS_STRING="F_MISSING>${MISSINGNESS}" #bcftools does not allow you to place a variable inside a string. The solution is to place the variable in a pre-made string, and refering to that instead.
MAF_STRING="MAX(AD/DP)>=${MAF} & MIN(AD/DP)<=$(awk -v maf=${MAF} 'BEGIN {printf "%.2f", 1-maf}')" #Same problem as the line above.
MAF_FILE_STRING=$(echo ${MAF} | sed -e 's/\.//g') #Make a string that does not include a '.' to create clean file names.
MISSINGNESS_FILE_STRING=$(echo ${MISSINGNESS} | sed -e 's/\.//g') #Make a string that does not include a '.' to create clean file names.
INT_DIR=${OUT_DIR}/IntermediateFiles_${FILENAME}_${MIN_RD}_Missing${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING} #Use a whole lot of information to specify intermediate directories.

#Create intermediate directory
mkdir -p ${INT_DIR}

#Change the version 4.3 vcf produced by VarScan to a 4.2 vcf.
zcat ${VCF_IN} | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/g' > ${INT_DIR}/${FILENAME}_V42.vcf

#Extract the mean depth across all pools per site using vcftools
vcftools --vcf ${INT_DIR}/${FILENAME}_V42.vcf --site-mean-depth --out ${INT_DIR}/ReadDepths

#Make sure that site depths with -nan read depth are ignored.
grep -v 'nan' ${INT_DIR}/ReadDepths.ldepth.mean > ${INT_DIR}/ReadDepthsCorrected.ldepth.mean

#use awk to get the mean and standard deviation in mean site read depth.
MEAN=$(tail -n +2 ${INT_DIR}/ReadDepthsCorrected.ldepth.mean | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
SD=$(tail -n +2 ${INT_DIR}/ReadDepthsCorrected.ldepth.mean | awk '{sum+=$3; sumsq+=$3*$3}END{print sqrt(sumsq/NR - (sum/NR)^2)}')

#Use the above values to calculate a threshold for maximum read depth filtering
MAX_RD_FLOAT=$(awk -v awkMEAN="${MEAN}" -v awkSD="${SD}" ' BEGIN { THRESHOLD=awkMEAN+2*awkSD; print THRESHOLD } ')
#This is not quite enough. bcftools can only take integers, so we need to round this number.
MAX_RD_INT=$(echo ${MAX_RD_FLOAT} | awk '{print int($1+0.5)}')

#Report the values in an intermediate file
echo "The average read depth calculated using vcftools --site-mean-depth and the standard deviation are as follows:" > ${INT_DIR}/FilteringSummary.txt
echo "Mean = ${MEAN}, SD = ${SD}" >> ${INT_DIR}/FilteringSummary.txt
echo "Making the maximum read depth threshold (Mean+2xSD) = ${MAX_RD_FLOAT}, which is rounded to ${MAX_RD_INT}" >> ${INT_DIR}/FilteringSummary.txt
echo "The maximum per-site missingness we are filtering for is ${MISSINGNESS}." >> ${INT_DIR}/FilteringSummary.txt
echo "The minor allele frequency filtering threshold is set to ${MAF}." >> ${INT_DIR}/FilteringSummary.txt

#Now that we have all the values we want, let's start filtering. First up is minimum and maximum read depth. We filter for this using VCFtools
vcftools --vcf ${INT_DIR}/${FILENAME}_V42.vcf --minDP ${MIN_RD} --max-meanDP ${MAX_RD_INT} \
--recode --recode-INFO-all --out ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}

#Next up, we filter out sites with poor quality (which won't filter out much), and filter for missingness using bcftools
bcftools +setGT ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}.recode.vcf -- -t q -n . -i 'GQ<15' | \
bcftools filter -e "${MISSINGNESS_STRING}" > ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}.vcf

# Filter for MAF>=5% in at least one population in the data set
bcftools view ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Missingness${MISSINGNESS_FILE_STRING}.vcf \
-i "${MAF_STRING}" -m2 \
> ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}.vcf

#There are two possible masks that can be applied. Some are .bed files with three columns specifying regions. Some are .txt files with two columns specifying single basepair sites (usually SNP locations).
#First, there are some lines that can apply a .txt mask. Then some lines that can apply .bed masks. Switching between them will necessitate editing this script.

#Applying .txt masks:
#Use bcftools to specify which sites to exclude.
bcftools view -T ^${HD_MASK} ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}.vcf \
> ${OUT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}_HDplot.vcf

# Produce summary of missingness per sample
vcftools --vcf ${OUT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}_HDplot.vcf \
--missing-indv --out ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}_HDplot

#Applying .bed masks:
#Create a .vcf file that contains only a header, and let bedtools write to this file. Bedtools doesn't output a header (for some reason).
#grep '#' ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Missingness${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}.vcf \
#> ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Missingness${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}_HDplotMask.vcf

#Finally, apply a mask created using HDplot output, that removes sites with excess heterozygosity.
#bedtools subtract -a ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Missingness${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}.vcf -b ${HD_MASK} \
#>> ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Missingness${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}_HDplotMask.vcf

#After all commands have run, clean up some unnecessary files left in the intermediate directory, and gzip up the useful intermediate vcfs.
rm ${INT_DIR}/${FILENAME}_V42.vcf
gzip ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}.recode.vcf
gzip ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}.vcf
gzip ${INT_DIR}/${FILENAME}_MinDP${MIN_RD}_MaxMeanDP${MAX_RD_INT}_Miss${MISSINGNESS_FILE_STRING}_MAF${MAF_FILE_STRING}.vcf

#As I was writing this script, I got insanely distracted at some point and made this chicken. Enjoy.

#   _\/_
#  / *  \      ////
# <      \____/  /   buh-gok
#  \            /     
#   \__     ___/     _
#      \___/        / \
#        |         |   |
#       //          \_/
