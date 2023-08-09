#!/bin/bash

# create output directory
mkdir -p ../results/align

# declare variables
IND=${1}
OUTPUT=../results/align/${IND}

# Remove PCR and optical duplicates with samtools
# This is the recommended way to mimic Picard tools MarkDuplicates with samtools
# See https://www.htslib.org/algorithms/duplicate.html

# Be mindful to account for as many cores as each piped command requires
# Here we have 4 single-threaded samtools commands so we need 4 cores
# samtools flagstat is run afterwards (much sorter time) utilizing all 4 cores

# sort reads in the BAM according to their names so that pairs are placed one below the other
# fix read-mates so that they both have the same sets of attributesfor the subsequent preprocessing
# sort reads again based on genomic coordinates and
# remove duplicates (-r argument)
samtools collate -O -u ${OUTPUT}.raw.bam \
    | samtools fixmate -m -u - - \
    | samtools sort -u - \
    | samtools markdup -r -f ${OUTPUT}_stats_file.txt -S -d 2500 --mode s --include-fails - ${OUTPUT}.markdup.bam

# gather statistics
# -@ number of cores
samtools flagstat -@4 ${OUTPUT}.markdup.bam > ${OUTPUT}.markdup.flagstat