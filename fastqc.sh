#!/bin/bash

#This script takes the name of a directory that contains fastq files and outputs a new directory in results that contains the fastqc reports. 

#Taking input
DIRECTORY_TO_QC=$1 # The full path to the directory that contains the fastq files you would like to perform QC on. 

# create output directory
OUTPUT_DIRECTORY=../results/fastqc_$(basename(${DIRECTORY_TO_QC})) #Use the name of the directory that stores the fastq files to name the output directory.
mkdir -p ${OUTPUT_DIRECTORY} # Make sure that the output directory exists.

# run fastqc
# edit number of cores used with the "-t" option
fastqc -t 4 -o ${OUTPUT_DIRECTORY} ${DIRECTORY_TO_QC}/*fastq.gz
