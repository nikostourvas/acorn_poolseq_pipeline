#!/bin/bash

# Trimmomatic
bash trimmomatic_simple_loop.sh

# Map
#bwa index /home/tourvasn/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa
#samtools faidx /home/tourvasn/ngs_training/data/reference/Plomion_et_al_2018/Qrob_V2_2N.fa

parallel --verbose -j 2 \
	'bash parallel_bwa.sh {}' :::: inds

# Picard
parallel --verbose -j 2 \
	'bash picard.sh {}' :::: inds

# Index all bam files again
samtools index -@ 8 /home/tourvasn/ngs_training/results/align/Plomion/*_sort.rmd.bam

# Variant calling
bash variant_calling_samtools-VarScan.sh
