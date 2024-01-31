#!/bin/bash

# variables
SAMPLE=$1
INPUT_DIR=/home/geneticsShare/AcornSeqData/AdapterClipped_Batch2
OUT_DIR=/home/geneticsShare/results/fastp_dedup_trim

mkdir -p ${OUT_DIR}

fastp --in1 ${INPUT_DIR}/${SAMPLE}_R1.fastq.gz --in2 ${INPUT_DIR}/${SAMPLE}_R2.fastq.gz \
      --stdout \
      --disable_adapter_trimming \
      --dedup --dup_calc_accuracy 3 \
      --length_required 50 \
      --correction \
      --html --report_title "${OUT_DIR}/${SAMPLE}_fastp_dedup.html" \
      --thread 2 |
fastp --stdin \
      --cut_right --cut_right_window_size 4 --cut_right_mean 20 \
      --html --report_title "${OUT_DIR}/${SAMPLE}_fastp_trim.html" \
      --out1 ${OUT_DIR}/${SAMPLE}_1.trim.dedup.fastq.gz --out2 ${OUT_DIR}/${SAMPLE}_2.trim.dedup.fastq.gz



      
