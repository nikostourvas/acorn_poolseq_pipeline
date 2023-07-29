#!/bin/bash

# Nikos Tourvas
# Convert FASTQ files compressed as bz2 to gz. The original files are retained.
# The script expects that each FASTQ file is found inside a separate directory.

INPUT=../data/untrimmed_fastq_ind/AdapterClipped
PATTERN=s*/*.bz2

cd $INPUT
find $PATTERN -type f > indfile

parallel --dry-run --verbose -j 18 \
    'bzcat {} | gzip -c > {.}.gz' :::: indfile

rm indfile
