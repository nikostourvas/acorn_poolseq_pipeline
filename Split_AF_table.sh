#!/bin/bash

#Lars Littmann
#20.06.2024
#Take an allele frequency table and divide it up into 30 equally sized, numbered files.
#This is done to parallelise lfmm.

AF_TABLE=$1

HEADER=$(head -n 1 ${AF_TABLE})

tail -n +2 ${AF_TABLE} > body.txt

split -n l/30 --numeric-suffixes=01 --additional-suffix .txt body.txt ${AF_TABLE/.txt/_Chunk}

sed -i "1i ${HEADER}" ${AF_TABLE/.txt/_Chunk??.txt}
