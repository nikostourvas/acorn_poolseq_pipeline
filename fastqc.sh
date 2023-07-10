#!/bin/bash

# Variables
READTYPE=$1 # either trimmed_fastq or untrimmed_fastq
DIRECTORY=$2 # either /mnt/data or ..results
# create output directory 
mkdir -p ../results/fastqc_$1

fastqc -t 32 -o ../results/fastqc_$1/ $2/$1/*fastq.gz
