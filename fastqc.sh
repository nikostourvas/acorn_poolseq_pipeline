#!/bin/bash

# Variables
READTYPE=$1 # either trimmed_fastq or untrimmed_fastq
DIRECTORY=$2 # either /mnt/data or ..results
# create output directory 
mkdir -p ../results/fastqc_$1

# run fastqc
# edit number of cores used with the "-t" option
fastqc -t 4 -o ../results/fastqc_$1/ $2/$1/*fastq.gz
