#!/bin/bash

# declare variables
IND=${1}
OUTPUT=../results/align_Batch1/${IND}

# index
samtools index ${OUTPUT}.markdup.Q20.bam