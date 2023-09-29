#!/bin/bash

#Lars Littmann
#29-09-2023
#This script creates a mappability mask using Heng Li's SNPable approach (Li, 2009). This script is adapted from a script provided by Axel Jensen, Uppsala University, department of ecology and genetics.
#The script relies on seqbility being installed (available at https://lh3lh3.users.sourceforge.net/snpable.shtml) and msmc-tools (available at https://github.com/stschiff/msmc-tools/tree/master)
#The only required input is a reference genome. Another thing to consider is the read length. (150 here, but may need to be changed depending on the application)

# refgenome as first input
GENOME=../reference/Qrob_PM1N.fa
# outdir as second
OUTDIR=../reference/SNPability_Mask/
#Temporary directory for intermediate files
INTDIR=../reference/SNPability_Mask/Intermediates
# set the path to where the seqbility scripts are 
SEQBILITY_PATH=../reference/SNPability_Mask/seqbility-20091110/
# tools to msmc-tools (used in the last step)
MSMC_TOOLS_PATH=../reference/SNPability_Mask/msmc-tools/

# create outdir if it doesn't exist
mkdir -p ${OUTDIR}

#Create directory for intermediate files
mkdir -p ${INTDIR}

# split the refgenome in to kmers as artificial reads, using 150 since that's our readlenght
${SEQBILITY_PATH}/seqbility-20091110/splitfa ${GENOME} 150 > ${INTDIR}/$(basename ${GENOME})_150splits.fa
# map them back to the reference chromosome
bwa aln -t 20 -R 1000000 -O 3 -E 3 ${GENOME} ${INTDIR}/$(basename ${GENOME})_150splits.fa > ${INTDIR}/$(basename ${GENOME})_150splits.aln.sai
# convert sai to sam
bwa samse ${GENOME} ${INT_DIR}/$(basename ${GENOME})_150splits.aln.sai ${INTDIR}/$(basename ${GENOME})_150splits.fa > ${INTDIR}/$(basename ${GENOME})_150splits.aln.sam
# compress it
bgzip ${INTDIR}/$(basename ${GENOME})_150splits.aln.sam
# generate the raw mask file
gzip -dc ${INTDIR}/$(basename ${GENOME})_150splits.aln.sam.gz | seqbility/seqbility-20091110/gen_raw_mask.pl > ${INTDIR}/rawMask_$(basename ${GENOME})_150.fa
# convert step 1
mkdir -p ${OUTDIR}/mappability_masks
${SEQBILITY_PATH}/seqbility-20091110/gen_mask -l 150 -r 0.5 ${INTDIR}/rawMask_$(basename ${GENOME})_150.fa > ${OUTDIR}/mappability_masks/mask_$(basename ${GENOME})_150.fa
# convert step 2, using the python script from msmc tools. Modified this to take input and output as command line arguments, for easier automation
# IMPORTANT: this script was written in python 2x something, it will not work under python3!

python2 ${MSMC_TOOLS_PATH}/makeMappabilityMask.py -i ${INTDIR}/mask_$(basename ${GENOME})_150.fa -o ${OUTDIR}/mappability_masks/mask_$(basename ${GENOME})_150.bed.gz &&

rm -r ${INTDIR}
