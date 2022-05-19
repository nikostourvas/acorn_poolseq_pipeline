#!/bin/bash

BASE="/home/tourvasn/ngs_training/results/trimmed_fastq"
RESULTS="flash_res"

# create directory for results
mkdir $RESULTS

# declare variables
IND=$1
FORWARD=$BASE/${IND}_R_1.trim.fastq.gz
REVERSE=$BASE/${IND}_R_2.trim.fastq.gz
OUTPUT=$RESULTS/${IND}.flash.merged

# run flash
flash -m 20 -x 0.05 -M 120 -o $OUTPUT -z -t 8 $FORWARD $REVERSE \
> $OUTPUT.flash.merged.out \
2 > $OUTPUT.flash.merged.err

