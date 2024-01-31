#!/bin/bash

# variables
SAMPLE=$1
INPUT_DIR=/home/geneticsShare/AcornSeqData/AdapterClipped_Batch2
OUT_DIR=/home/geneticsShare/AcornSeqData/trimmed

mkdir -p ${OUT_DIR}

fastp --in1 ${INPUT_DIR}/${1}_R1.fastq.gz --in2 ${INPUT_DIR}/${1}_R2.fastq.gz \
      --out1 ${OUT_DIR}/${1}_1.trim.dedup.fastq.gz --out2 ${OUT_DIR}/${1}_2.trim.dedup.fastq.gz \
      --disable_adapter_trimming \
      --dedup --dup_calc_accuracy 3 \
      --length_required 50 \
      --correction \
      --html --report_title "${OUT_DIR}/${1}_fastp.html" \
      --thread 2