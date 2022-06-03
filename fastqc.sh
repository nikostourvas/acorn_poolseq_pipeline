#!/bin/bash

# Variables
READTYPE=$1 # either trimmed_fastq or untrimmed_fastq
DIRECTORY=$2 # either data or results
# create output directory 
mkdir ../results/fastqc_$1

fastqc -t 4 -o ../results/fastqc_$1/ ../$2/$1/*fastq.gz
