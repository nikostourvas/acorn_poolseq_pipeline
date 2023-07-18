#!/bin/bash

#Lars Littmann
#26.06.2023
#Comparing allele frequencies of a population determined from pool sequencing data with allele frequencies determined from 
#sequences of individuals of a population by linear regression of the allele frequencies. Dependent on R script. See Github

#INPUTS:
#1.Hard-filtered VCF created from pool-sequencing data. No need for read depth filtering, the script does this. 
#Unpruned data. The relevant pool needs to be specified by modifying this script.
#2.Hard-filtered VCF created from individual-sequencing data. Should contain at least 20 individuals that can be compared against one pool. 
#No need for Read depth filtering or pruning.
#3.A numeric value for the threshold of minimum read depth per individual
#4.A numeric value for the threshold of minimum cumulative read depth of the pool

#PROCESSING VCF FILES BEFORE STARTING IN R

#LINING UP ALL OF OUR INPUTS.

#When running this script from the command line, specify the following three inputs:
POOL_VCF=$1 #(Gzipped) VCF file of PoolSeq data. I used home/llittmann/local_disks/disk_3/LarsAlleleFreqVal/allRegs_Qpe_all_targets_intervals_AquTar_sorted.vcf.gz
IND_VCF=$2 #(Gzipped) VCF file of IndSeq data
IND_MIN_RD=$3 # The minimum readdepth that EACH individual should have to be part of the IndSeq data.
R_SCRIPT=./AlleleFreqRegression_SupportingScript.R # The R script that will be used to further filter allele frequencies and perform the regression.
#The way I direct the script to the location of this R script is very sloppy, but since everyone is working on a seperate container it's the best I can do.

#CREATING THE NECESSARY DIRECTORIES.

#Create a subdirectory of the current working directory and name it using the filenames of the VCFs.
OUTDIR=./AlleleFreqRegression_$(basename ${POOL_VCF/.vcf.gz/})_WITH_$(basename ${IND_VCF/.vcf.gz/})
mkdir -p ${OUTDIR} #Create the directory if it does not exist yet.
#Also create a directory within the working directory to place intermediate nonsense.
INTERMEDIATEDIR=/tmp

#The following two lines can be used instead of '/tmp', in case containers do not support this functionality.
#./AlleleFreqRegression_$(basename ${POOL_VCF/.vcf.gz/})_WITH_$(basename ${IND_VCF/.vcf.gz/})/IntermediateFiles
#mkdir -p ${INTERMEDIATEDIR}


#PROCESSING VCF FILES BEFORE IMPORTING IN R

#For PoolSeq data

#Following one-liner extracts allele frequency for each site on the first pool. This should yield all the info necessary to run the R-script later.
bcftools query -s Sample4 -f '%CHROM %POS %REF %ALT [%DP %AD %FREQ]\n' ${POOL_VCF} > ${INTERMEDIATEDIR}/Pool_Sample1_AltCalls_Freq.txt

#For IndSeq data

#The following vcftools command first excludes GENOTYPES that are not based on the minimum read depth (--minDP).
#It then excludes SITES for which one or more genotypes were removed in the last step (--max-missing-count 0).
#Not sure why, but I also encountered SNPs that were not biallelic, so I also enforce that by setting the minimum and maximum  number of alleles to 2 (--min-alleles and --max-alleles)
vcftools --gzvcf ${IND_VCF} --minDP ${IND_MIN_RD} --max-missing-count 0 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > ${INTERMEDIATEDIR}/$(basename ${IND_VCF/.vcf.gz/_RDfiltered.vcf.gz});

#Next, use vcftools again on the output data and create a file with allele frequencies. Seems ideal, but the format it comes in is quite annoying.
vcftools --gzvcf ${INTERMEDIATEDIR}/$(basename ${IND_VCF/.vcf.gz/_RDfiltered.vcf.gz}) --freq --out ${INTERMEDIATEDIR}/$(basename ${IND_VCF/.vcf.gz/_AlleleFreq});

#This .freq file has a header that messes with R, so we remove it using tail and pipe it into a file with a better name.
tail -n +2 ${INTERMEDIATEDIR}/$(basename ${IND_VCF/.vcf.gz/_AlleleFreq.frq}) > ${INTERMEDIATEDIR}/Processed_IndSeq_ForR.txt

#PASSING ON TO R

#This line passes on the two processed VCFs (Pool and Ind) to the R-script that further filters, merges the data, performs regression, and outputs figures.
Rscript ${R_SCRIPT} ${INTERMEDIATEDIR}/Pool_Sample1_AltCalls_Freq.txt ${INTERMEDIATEDIR}/Processed_IndSeq_ForR.txt ${OUTDIR} #Rscript runs the R file and passes on the two VCF files as arguments. Third argument is the out directory

#TESTING THE SCRIPT